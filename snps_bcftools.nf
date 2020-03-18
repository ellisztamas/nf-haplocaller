/*
BCFtools to call SNPs in the known sites.


nextflow run snpsnf --input "*bam" --targets target_file --fasta ref_seq/TAIR10_wholeGenome.fasta  --outdir filter_vcfs
*/

params.project = "the1001genomes"
build_index = false
params.outdir = './snpcall'
params.fasta = false

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

input_bams = Channel
      .fromPath( "${params.input}" )
      .map { it -> [ file("$it").baseName, file("$it"), file("${it}.bai") ] }

target_file = Channel.fromPath( "$params.targets" )

all_inputs = input_bams.combine(target_file)


process getVcf {
  tag "$name"
  label 'env_medium'
  publishDir "$params.outdir", mode: "copy"

  input:
  set val(name), file(bam), file(bam_index), file(known_sites) from all_inputs
  // set file(tar_file), file(tar_index) from target_file.collect()

  output:
  file("${name}.vcf*") into output_vcf

  script:
  // .sorted.mkdup
  """
  bcftools mpileup -R $known_sites -f $reffol/${refid}.fasta $bam | bcftools call -c -T $known_sites | bcftools filter -e 'QUAL < 20' -O z -o ${name}.vcf.gz
  """
  // #bcftools view -T $known_sites $bam -O v | bgzip -c > ${name}.vcf.gz
  // #bcftools index ${name}.vcf.gz
  // bcftools mpileup -Ou -f $reffol/${refid}.fasta $bam | bcftools call -Ov -T $known_sites -C alleles -m  > ${name}.vcf
}
