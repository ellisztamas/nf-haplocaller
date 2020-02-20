/*
Provide directories as inputs with .gvcf files present

nextflow run ~/mygit/nf-haplocaller/get_gvcf.nf --input "0*"  --input_ext "gvcf" --fasta 100.REF_SEQ/TAIR10_wholeGenome.fasta --outdir temp

*/

params.input = false
params.input_ext = "gz"

params.project = "the1001genomes"
build_index = false
params.outdir = '.'
params.fasta = false

params.cohort = "allsample"
params.name = false
chromosomes = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]

if ( params.fasta ){
  genome = file(params.fasta)
  reffol = genome.parent
  refid = genome.baseName
  if( !genome.exists() ) exit 1, "Reference fasta file not found: ${params.fasta}"
  bwa_indices = Channel
    .fromPath( "$reffol/${refid}.fasta.b*" )
    .ifEmpty { build_index = true }
} else {
  exit 1, "Provide reference fasta file. Ex., --fasta file_path"
}

/*
* 1. Create a channel for checking bwa index for genome ref
*/
if (build_index == true){
  process makeBWAindex {
      publishDir "${reffol}", mode: 'copy'

      input:
      file genome

      output:
      set  file("${refid}.fasta"), file("${refid}.fasta.*"), file "${refid}.dict" into fasta_index

      script:
      """
      samtools faidx ${genome}
      bwa index $genome
      picard  CreateSequenceDictionary -R $genome -O ${refid}.dict
      """
  }
} else {
  fasta_index = Channel
    .from( [file("$reffol/${refid}.fasta"), file("$reffol/${refid}.fasta.*"), file("${refid}.dict") ])
}

/*
* 10.1 GenomicsDBImport for all the files
*/
chromosomes_ch = Channel.from( chromosomes )

raw_gvcf =  Channel.fromPath( "${params.input}" )
raw_gvcf_index = Channel.fromPath( "${params.input}.*" )


process GenomicsDBImport {
  tag "gendbi_${chr}"
  label 'env_large'
  publishDir "$params.outdir/raw_variants/genodbi_gatk", mode: 'copy'

  input:
  each chr from chromosomes_ch
  file in_vcf from raw_gvcf.collect()
  file(vcf_idx) from raw_gvcf_index.collect()

  output:
  set val(chr), file ("${params.cohort}.genomicsdbi.${chr}.dbi") into gendbi_in

  script:
  def try_vcfs = in_vcf.collect { "-V $it" }.join(' ')
  """
  gatk GenomicsDBImport --java-options "-Xmx${task.memory.toGiga()}G -Xms${task.memory.toGiga()}G" \
  --tmp-dir=${params.tmpdir} \
  ${try_vcfs} \
  -L ${chr}\
  --batch-size 50 \
  --genomicsdb-workspace-path ${params.cohort}.genomicsdbi.${chr}.dbi
  """
}

/*
* 10.2 GenotypeGVCF for all the files
*/
process genotypeGVCFs {
  tag "gvcf_${chr}"
  label 'env_large'
  // publishDir "$params.outdir/combinedGVCF", mode: 'copy'

  input:
  set val(chr), file(workspace) from gendbi_in

  output:
  set val(chr), file("${params.cohort}.${chr}.vcf.gz"), file("${params.cohort}.${chr}.vcf.gz.tbi") into combined_vcf

  script:
  """
  WORKSPACE=\$( basename ${workspace} )
  gatk GenotypeGVCFs --java-options "-Xmx24g -Xms24g" \
  --tmp-dir ${params.tmpdir}\
  -R $reffol/${refid}.fasta \
  -O ${params.cohort}.${chr}.vcf.gz \
  -V gendb://\$WORKSPACE \
  -L $chr \
  --only-output-calls-starting-in-intervals \
  --use-new-qual-calculator
  """
  // -D ${dbsnp_resource_vcf} \
  // -G StandardAnnotation \
}

/*
* 11 Join all the VCFs generated per chromosome
*/
try_vcf = chromosomes.collect{it -> "--INPUT ${params.cohort}.${it}.vcf.gz"}.join(' ')
process gatherVcfs {
  tag "gatherVCF_${params.cohort}"
  label 'env_gatk_medium'
  publishDir "$params.outdir/raw_variants/combined_GVCF", mode: 'copy'

  input:
  file (vcf) from combined_vcf.collect()

  output:
  set file("${params.cohort}.vcf.gz"), file("${params.cohort}.vcf.gz.tbi") into (gvcf_sample_ch, gvcf_select_snp_ch)

  script:
  """
  gatk --java-options "-Xmx10g -Xms10g" \
  GatherVcfs \
  ${try_vcf} \
  --OUTPUT ${params.cohort}.vcf.gz
  tabix ${params.cohort}.vcf.gz
  """
}

/*
* 12.1 Get sample names and filter the GVCF by SelectVariants
*/
process getSamples {
  tag "samples_${params.cohort}"
  label 'env_small'

  input:
  set file(gvcf), file(gvcf_index) from gvcf_sample_ch

  output:
  stdout sample_names

  script:
  """
  tabix -H $gvcf | grep -m 1 "^#CHROM" | awk '{for(i=10;i<=NF;i++) print \$i}'
  """
}

input_names = sample_names
    .splitCsv()
    .map {row -> "${row[0]}"}


process selectSNPs {
  tag "$name"
  publishDir "${params.outdir}/raw_variants/sample_biallelics", mode: 'copy'
  label 'env_gatk_medium'

  input:
  val name from input_names
  set file(gvcf), file(gvcf_index) from gvcf_select_snp_ch.collect()

  output:
  file("${name}.BIALLELIC.SNPs.vcf.gz*") into filter_vcf

  script:
  """
  gatk SelectVariants --java-options "-Xmx24g -Xms24g" \
    --tmp-dir ${params.tmpdir}\
    -R $reffol/${refid}.fasta \
    -V $gvcf \
    -select-type SNP --restrict-alleles-to BIALLELIC \
    -se "${name}" \
    -O ${name}.BIALLELIC.SNPs.vcf.gz
  """
}



// input_gvcfs = Channel
//       .fromPath( "${params.input}", type: 'dir' )
//       .map{ it -> [it.name, file("$it/*${params.input_ext}"), file("$it/*${params.input_ext}.*")] }
//
// process joinGVCFs {
//   label 'env_large'
//   tag "$fol_name"
//   // publishDir "$params.outdir", mode: 'copy'
//
//   input:
//   set val(fol_name), file(in_vcf), file(in_vcf_idx) from input_gvcfs
//   set ref, ref_idx, ref_dict from fasta_index.collect()
//
//   output:
//   set val("$fol_name"), file("${fol_name}.vcf.gz"), file("${fol_name}.vcf.gz.tbi") into combgVCF
//
//   script:
//   def try_vcfs = in_vcf.collect { "-V $it" }.join(' ')
//   """
//   java -Djava.io.tmpdir=${params.tmpdir} -jar \$EBROOTGATK/GenomeAnalysisTK.jar\
//   -T GenotypeGVCFs -R $ref\
//   -nt ${task.cpus} \
//   $try_vcfs -o ${fol_name}.vcf.gz \
//   -allSites
//   """
//   //gatk --java-options '-Djava.io.tmpdir=${params.tmpdir}' GenotypeGVCFs\
//   //-R $reffol/${refid}.fasta\
//   //$try_vcfs -O ${fol_name}.vcf.gz \
//   //- -variant_index_type LINEAR -variant_index_parameter 128000
// // }
//
// process selectSNPs {
//   tag "$fol_name"
//   label 'env_medium'
//   publishDir "$params.outdir", mode: 'copy'
//
//   input:
//   set val(fol_name), file(gvcf), file(gvcf_index) from combgVCF
//
//   output:
//   set file("${fol_name}.filter.vcf"), file("${fol_name}.filter.vcf.idx") into filter_vcf
//
//   script:
//   """
//   java -Djava.io.tmpdir=${params.tmpdir} -jar \$EBROOTGATK/GenomeAnalysisTK.jar\
//   -T SelectVariants -R $reffol/${refid}.fasta\
//   -nt ${task.cpus} \
//   -V $gvcf \
//   -o ${fol_name}.filter.vcf \
//   -selectType SNP -restrictAllelesTo BIALLELIC
//   """
// }
