/*
BCFtools to call SNPs in the known sites.


nextflow run snpsnf --input "*bam" --targets target_file --fasta ref_seq/TAIR10_wholeGenome.fasta  --outdir filter_vcfs
*/

params.project = "the1001genomes"
params.input = false
params.outdir = './snpcall'
params.fasta = false
// to perform snpcalling only on specific position, provide a VCF file below
params.known_sites = false  // if known sites (VCF file) is given,
// Provide the --cohort, if you want to generate a combined vcf for all the samples.
params.cohort = false //

// Options to filter VCF file
params.minimum_depth = 5
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

req_label = params.cohort != false ? 'env_large' : 'env_medium'
final_outdir = params.cohort != false ? params.outdir : "${params.outdir}/perSample"

process getVcf {
  tag "${outname}"
  label "${req_label}"
  publishDir "${final_outdir}", mode: "copy"

  input:
  file(bam) from input_bams
  file(bam_index) from input_bam_idx

  output:
  set val("${outname}"), file("${outname}.vcf*") into nofilter_vcf

  script:
  known_sites_cmd = params.known_sites != false ? "-T ${known_sites_vcf}" : ''
  outname = params.cohort != false ? "${params.cohort}" : file(bam).baseName
  """
  bcftools mpileup --threads ${task.cpus} $known_sites_cmd -f $reffol/${refid}.fasta $bam | bcftools call --threads ${task.cpus} -c $known_sites_cmd -O z -o ${outname}.vcf.gz
  """
}

process filterVcf {
  tag "${name}"
  label 'env_medium'
  publishDir "$params.outdir", mode: "copy"

  input:
  set val(name), file(vcf) from nofilter_vcf

  output:
  file("${name}.qualfiltered.vcf.gz") into output_vcf

  script:
  // qual_filter = $params.minimum_qual > 0 ? "  "
  """
  bcftools filter -e 'QUAL < $params.minimum_qual' -O z -o ${name}.qualfiltered.vcf.gz $vcf
  """
}
