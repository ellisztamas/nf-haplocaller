/*
 to call SNPs using bcftools either in known sites (cohort or not).
 you can also filter snps in case of bsseq

nextflow run snpsnf --input "*bam" --bsseq --cohort f2_snps --known_sites target_file --fasta ref_seq/TAIR10_wholeGenome.fasta  --outdir filter_vcfs
*/

params.project = "the1001genomes"
params.input = false
params.outdir = './snpcall'
params.fasta = false
// to perform snpcalling only on specific position, provide a VCF file below
params.known_sites = false  // if known sites (VCF file) is given,
// Provide the --cohort, if you want to generate a combined vcf for all the samples.
params.cohort = false // file name for output combined vcf
params.call_variants = true // Only call variants instead of all sites, useful if given cohort and known_sites (to then remove non-polymorphic in the sites)
// Can we filter SNPs for bisulfite
params.bsseq = false 
params.minimum_depth = 3

// Options to filter VCF file
params.minimum_qual = 30

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


//  SNP calling starts now
if (params.known_sites != false){
  known_sites_vcf = file( params.known_sites )
}

if (params.cohort != false) {
  input_bams = Channel.fromPath( "${params.input}" )
  input_bam_idx = Channel.fromPath( "${params.input}" ).map{ it -> file("${it}.bai")  }.collect()
  input_bams = input_bams.collect()
} else {
  input_bams = Channel.fromPath( "${params.input}" )
  input_bam_idx = Channel.fromPath( "${params.input}" ).map{ it -> file("${it}.bai")  }
}

call_varonly = params.call_variants != false ? '-v' : ' '
req_label = params.cohort != false ? 'env_large' : 'env_medium'

process getVcf {
  tag "${outname}"
  label "${req_label}"
  publishDir "${params.outdir}/raw_variants", mode: "copy"

  input:
  file(bam) from input_bams
  file(bam_index) from input_bam_idx

  output:
  set val("${outname}"), file("${outname}.vcf.gz") into nofilter_vcf

  script:
  known_sites_cmd = params.known_sites != false ? "-T ${known_sites_vcf}" : ''
  outname = params.cohort != false ? "${params.cohort}" : file(bam).baseName
  """
  bcftools mpileup --threads ${task.cpus} $known_sites_cmd -f $reffol/${refid}.fasta $bam | bcftools call --threads ${task.cpus} -c $call_varonly $known_sites_cmd -O z -o ${outname}.vcf.gz
  """
}

process filterVcf {
  tag "${name}"
  label 'env_medium'
  publishDir "${params.outdir}/filtered", mode: "copy"

  input:
  set val(name), file(vcf) from nofilter_vcf

  output:
  set val(name), file("${name}.filtered.vcf.gz"), file("${name}.filtered.vcf.gz.tbi") into output_vcf

  script:
  // qual_filter = $params.minimum_qual > 0 ? "  "
  """
  bcftools filter -e 'QUAL < $params.minimum_qual' -O z -o ${name}.filtered.vcf.gz $vcf
  tabix ${name}.filtered.vcf.gz
  """
}

if (params.bsseq != false){
  process bsseqVcf {
    tag "${name}"
    label 'env_medium'
    publishDir "${params.outdir}/bsseq", mode: "copy"

    input:
    set val(name), file(vcf), file(vcf_idx) from output_vcf

    output:
    set val(name), file("${name}.bsseq.filtered.bed") into bsseq_vcf

    script:
    // qual_filter = $params.minimum_qual > 0 ? "  "
    // echo "CHROM,POS,REF,ALT,DP," > ${name}.bsseq.filtered.bed
    // bcftools query -l $vcf | tr '\n' ',' | sed 's/,$/\n/' >> ${name}.bsseq.filtered.bed
    """
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%DP[\t%GT]\n" $vcf | awk '\$3 !~ /C|G/ && length(\$3) == 1 && length(\$4) == 1 && \$4 !~ /T/ &&  \$5 > $params.minimum_depth {print \$0 }' > ${name}.bsseq.filtered.bed
    """
  }
}

// process getSamples {
//   tag "joinGVCF"
//   label 'env_gatk_small'
//
//   input:
//   file gvcf from combgVCF_name
//   file gvcf_index from combgVCF_name_index
//
//   output:
//   stdout sample_names
//
//   script:
//   """
//   tabix -H $gvcf | grep -m 1 "^#CHROM" | awk '{for(i=10;i<=NF;i++) print \$i}'
//   """
// }
//
// input_names = sample_names
//     .splitCsv()
//     .map {row -> "${row[0]}"}
//
