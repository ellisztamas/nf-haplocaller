/*
 to call SNPs using BisSNP in known sites 

nextflow run snpsnf --input "*bam" --known_sites 1001g_vcf --fasta ref_seq/TAIR10_wholeGenome.fasta  --outdir filter_vcfs
*/

params.project = "the1001genomes"
params.input = false
params.outdir = './snpcall'
params.fasta = false
// to perform snpcalling only on specific position, provide a targets file
params.known_sites_vcf = "/groups/nordborg/projects/the1001genomes/scratch/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.vcf"
// known sites (VCF file) is given,
params.known_sites_targets = "/groups/nordborg/projects/the1001genomes/scratch/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.tsv.gz"  
// targets file generated from bcftools

// build_index = false
// if ( params.fasta ){
//   genome = file(params.fasta)
//   reffol = genome.parent
//   refid = genome.baseName
//   if( !genome.exists() ) exit 1, "Reference fasta file not found: ${params.fasta}"
//   bwa_indices = Channel
//     .fromPath( "$reffol/${refid}.fasta.b*" )
//     .ifEmpty { build_index = true }
// } else {
//   exit 1, "Provide reference fasta file. Ex., --fasta file_path"
// }


//  SNP calling starts now
input_bams = Channel
        .fromPath( "${params.input}" )
        .map{ it -> [ file(it).baseName, file(it), file("${it}.bai") ] }
known_sites_vcf = file( params.known_sites_vcf  )
known_sites_targets = file( params.known_sites_targets  )


process ModifyBam {
  tag "${name}"
  label "env_medium"
  storeDir "${params.outdir}/modified_bams"

  input:
  set val(name), file(bam), file(bam_index) from input_bams

  output:
  set val(name), file("${name}.modified.bam") into modify_bam

  script:
  """
  picard AddOrReplaceReadGroups I=$bam \
  O=${name}.modified.bam \
  RGID=$name RGLB=$name \
  RGPL=illumina RGPU=none RGSM=$name
  """
}


process getBSVcf {
  tag "${name}"
  label "env_medium"
  publishDir "${params.outdir}/raw_variants", mode: "copy"

  input:
  set val(name), file(bam) from modify_bam

  output:
  file("${name}.filtered.snpvcf.bed") into bsseq_vcfbed

  script:
  outname = bam.baseName 
  """
  samtools index $bam 
  gatk --java-options '-Djava.io.tmpdir=${params.tmpdir}' HaplotypeCaller\
  -R $reffol/${refid}.fasta\
  -I $bam \
  -O ${name}.vcf.gz \
  --alleles ${known_sites_vcf} \
  --force-call-filtered-alleles \
  --output-mode EMIT_ALL_CONFIDENT_SITES

  bcftools filter -T $known_sites_targets ${name}.vcf.gz | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%AD][\t%GT]\n" |\
  awk 'length(\$3) == 1 && length(\$4) == 1 {print \$0}'  > ${name}.snpvcf.bed
  
  python $workflow.projectDir/scripts/01.filter_bsseq_variants.py \
  -i ${name}.snpvcf.bed -o ${name}.filtered.snpvcf.bed
  """
}

// process getBSVcf {
//   tag "${name}"
//   label "env_medium"
//   publishDir "${params.outdir}/raw_variants", mode: "copy"

//   input:
//   set val(name), file(bam) from modify_bam

//   output:
//   file("${name}.genotype.vcf.gz") into bsseq_vcf

//   script:
//   outname = bam.baseName
//   """
//   samtools index $bam 
//   java -Xmx20G -jar $workflow.projectDir/BisSNP/BisSNP-1.0.0.jar \
//   -R $reffol/${refid}.fasta \
//   -I $bam\
//   -D ${known_sites_vcf} \
//   -T BisulfiteGenotyper \
//   -vfn1 ${name}.genotype.raw.vcf -vfn2 ${name}.modified.snp.raw.vcf \
//   -C CG,1 -C CH,1 -out_modes EMIT_ALL_SITES \
//   -stand_call_conf 10 -minConv 1 -vcfCache 1000000 -mmq 30 -mbq 5

//   bcftools filter -T $known_sites_targets\
//   ${name}.genotype.raw.vcf \
//   -O z -o ${name}.genotype.vcf.gz
//   """
// }

    