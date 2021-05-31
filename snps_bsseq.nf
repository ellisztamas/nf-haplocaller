/*
 to call SNPs using BisSNP in known sites 

nextflow run snpsnf --input "*bam" --known_sites 1001g_vcf --fasta ref_seq/TAIR10_wholeGenome.fasta  --outdir filter_vcfs
*/

params.project = "the1001genomes"
params.input = false
params.outdir = './snpcall'
params.fasta = false
params.cohort = false // file name for output combined vcf
// add additional information in the VCF file. check 
// http://samtools.github.io/bcftools/bcftools.html#mpileup
params.annotate = "FORMAT/DP"

/* Below are additional filtering steps, only work on cohort analysis
1. filter sites based on minor allele frequency (maf)
2. missing fraction of missing calls
*/
params.filter_maf = false
params.min_maf = 0.2  //
params.max_na = 0.3 //

params.snpcaller = "bcftools"
// You could also give gatk, bcftools

params.min_base_quality = 30 // only used for snpcaller bcftools
params.min_snp_qual = 30
// to perform snpcalling only on specific position, provide a targets file
// params.known_sites_vcf = "/groups/nordborg/projects/the1001genomes/scratch/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.vcf"
// known sites (VCF file) is given,
params.known_sites = "/groups/nordborg/projects/the1001genomes/scratch/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.tsv.gz"  
// targets file generated from bcftools

build_index = false
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


ch_input_raw_bams = Channel
      .fromPath( "${params.input}" )
      .map{ it -> [ file(it).baseName, file(it) ] }
/*
//  SNP calling after modifying the bam file

Please cite: 
Adam Nunn et. al. 2021
https://www.biorxiv.org/content/10.1101/2021.01.11.425926v1.full.pdf

*/

process modifyBam {
  tag "${sample_id}"
  storeDir "${workflow.workDir}/modifiedBam"
  label "env_medium"

  input:
  set val(name), file(bam) from ch_input_raw_bams

  output:
  set val(sample_id), file("${sample_id}.bsmod.bam") into ch_modify_bam
  set val(sample_id), file("${sample_id}.bsmod.bam.bai") into ch_modify_bam_index

  script:
  sample_id = name.replaceAll( "_1_val_1_bismark_bt2_pe.deduplicated", "")
  sample_id = sample_id.replaceAll("_1.sorted.markDups", "")
  def avail_mem = task.memory ? ((task.memory.toGiga() - 6) / task.cpus).trunc() : false
  def sort_mem = avail_mem && avail_mem > 2 ? "-m ${avail_mem}G" : ''
  """
  samtools sort $sort_mem --threads ${task.cpus} -o ${sample_id}.sorted.bam $bam

  samtools calmd -b ${sample_id}.sorted.bam $reffol/${refid}.fasta 1> ${sample_id}_calmd.bam 2> samtools.log.err
  samtools index ${sample_id}_calmd.bam

  python $workflow.projectDir/scripts/epidiverse_change_sam_queries.py \
    -f $reffol/${refid}.fasta -T ${task.cpus} -Q \
    -t ${workflow.workDir}/ \
    ${sample_id}_calmd.bam  ${sample_id}.raw.bsmod.bam

  picard AddOrReplaceReadGroups\
  I=${sample_id}.raw.bsmod.bam O=${sample_id}.bsmod.bam\
  ID=$sample_id LB=$sample_id PL=illumina PU=none SM=$sample_id

  samtools index ${sample_id}.bsmod.bam
  """
}

/*
Using BCFtools to call genotypes at known sites
*/

if (params.snpcaller == "bcftools"){

  process BCFcall {
    tag "${name}"
    label "env_medium"
    publishDir "${params.outdir}/variants_bcftools", mode: "copy"

    input:
    set val(name), file(bam) from ch_modify_bam
    set val(name), file(bam_index) from ch_modify_bam_index

    output:
    file("${name}.snps.vcf.gz*") into ch_bsseq_snps

    script:
    outname = bam.baseName 
    // flagW=99,147 flagC=83,163 // set for paired end
    known_sites_mpile_cmd = params.known_sites != false ? "-T ${params.known_sites}" : ''
    known_sites_call_cmd = params.known_sites != false ? "-C alleles -m -T ${params.known_sites}" : '-m'
    """
    bcftools mpileup --threads ${task.cpus}\
    $known_sites_mpile_cmd \
    -a $params.annotate \
    -Q ${params.min_base_quality}\
    -f $reffol/${refid}.fasta $bam | \
    bcftools call --threads ${task.cpus} \
    $known_sites_call_cmd \
    -O z -o ${name}.raw.vcf.gz

    bcftools filter -e 'QUAL < $params.min_snp_qual' \
    -Oz ${name}.raw.vcf.gz | \
    bcftools view -V indels > ${name}.snps.vcf
    bgzip ${name}.snps.vcf && tabix ${name}.snps.vcf.gz
    
    """
    // bcftools filter -e ' REF == "G" && ALT == "A" ' ${name}.snps.vcf.gz |\
    // bcftools filter -e ' REF == "C" && ALT == "T" ' > ${name}.qual_filtered.vcf
    // bgzip ${name}.qual_filtered.vcf && tabix ${name}.qual_filtered.vcf.gz
  }

}

/*
Using GATK haplotyper caller below to call SNPs
*/

if (params.snpcaller == "gatk"){
  
  process gatkCaller {
    // Using GATK
    tag "${name}"
    label "env_medium"
    publishDir "${params.outdir}/variants", mode: "copy"

    input:
    set val(name), file(bam) from ch_modify_bam
    set val(name), file(bam_idx) from ch_modify_bam_index

    output:
    file("${name}.vcf.gz*") into ch_bsseq_snps

    script:
    known_sites_cmd = params.known_sites != false ? "--alleles ${params.known_sites} --force-call-filtered-alleles" : ''
    """
    gatk --java-options '-Djava.io.tmpdir=${params.tmpdir}' HaplotypeCaller\
    -R $reffol/${refid}.fasta\
    -I $bam \
    -O ${name}.raw.vcf.gz \
    $known_sites_cmd \
    --output-mode EMIT_ALL_CONFIDENT_SITES

    bcftools filter -e 'FILTER == "LowQual"' -O z ${name}.raw.vcf.gz |\
    bcftools view -V indels > ${name}.vcf
    bgzip ${name}.vcf && tabix ${name}.vcf.gz
    """
  }
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

if (params.cohort != false) {

  process mergeVCF {
    label "env_large"
    publishDir "${params.outdir}/", mode: "copy"

    input:
    file(vcf) from ch_bsseq_snps.collect()

    output:
    file("${params.cohort}.qual_filtered.vcf.gz") into mergedVCF

    script:
    """
    bcftools merge *.vcf.gz -Oz -o ${params.cohort}.qual_filtered.vcf.gz
    """
  }
  if (params.filter_maf) {
    process filterMAF {
      label "env_medium"
      publishDir "${params.outdir}/", mode: "copy"

      input:
      file(vcf) from mergedVCF

      output:
      file("${params.cohort}.maf_filtered.vcf.gz*") into filter_maf_vcf

      script:
      """
      bcftools view -i 'F_MISSING <= $params.max_na && MAF>=$params.min_maf' $vcf > ${params.cohort}.maf_filtered.vcf
      bgzip ${params.cohort}.maf_filtered.vcf && tabix ${params.cohort}.maf_filtered.vcf.gz 
      """
    }
  }
  
} 
