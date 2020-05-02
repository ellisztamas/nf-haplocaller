/*
 to call SNPs using bcftools either in known sites (cohort or not).

nextflow run snpsnf --input "*bam" --variants_only --cohort f2_snps --known_sites target_file --fasta ref_seq/TAIR10_wholeGenome.fasta  --outdir filter_vcfs
*/

params.project = "the1001genomes"
params.input = false
params.outdir = './snpcall'
params.fasta = false
// provide save_intermediate if you want to save vcf file without filtering
params.save_intermediate = false
// to perform snpcalling only on specific position, provide a targets file
params.known_sites = false  // if known sites (VCF file) is given,
// Provide the --cohort, if you want to generate a combined vcf for all the samples.
params.cohort = false // file name for output combined vcf
params.variants_only = false // Only call variants instead of all sites, useful if given cohort and known_sites (to then remove non-polymorphic in the sites)
// Can we filter SNPs for bisulfite
params.get_bed = false  // will output only bed file 

// Below options important to consider while calling SNPs in different samples
params.min_depth = 2
params.min_base_quality = 30
params.min_map_quality = 10 // really important > 10 would be for uniquly mapped reads
//http://seqanswers.com/forums/showthread.php?t=61908
// Options to filter VCF file
params.min_snp_qual = 30


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
  known_sites_tsv = file( params.known_sites )
  //     bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $vcf |\
  //     bgzip -c > ${out_name}.tsv.gz && \
  //     tabix -s1 -b2 -e2 ${out_name}.tsv.gz
}

if (params.cohort != false) {
  input_bams = Channel.fromPath( "${params.input}" )
  input_bam_idx = Channel.fromPath( "${params.input}" ).map{ it -> file("${it}.bai")  }.collect()
  input_bams = input_bams.collect()
} else {
  input_bams = Channel.fromPath( "${params.input}" )
  input_bam_idx = Channel.fromPath( "${params.input}" ).map{ it -> file("${it}.bai")  }
}

call_varonly = params.variants_only != false ? '-v' : ' '
req_label = params.cohort != false ? 'env_large' : 'env_medium_mem'

process getVcf {
  tag "${outname}"
  label "${req_label}"
  publishDir "${params.outdir}/raw_variants", mode: "copy",
    saveAs: { filename -> params.save_intermediate ? filename : null }

  input:
  file(bam) from input_bams
  file(bam_index) from input_bam_idx

  output:
  set val("${outname}"), file("${outname}.vcf.gz") into nofilter_vcf

  script:
  known_sites_mpile_cmd = params.known_sites != false ? "-T ${known_sites_tsv}" : ''
  known_sites_call_cmd = params.known_sites != false ? "-C alleles -m -T ${known_sites_tsv}" : ''
  outname = params.cohort != false ? "${params.cohort}" : file(bam).baseName
  // -q ${params.min_map_quality}\
  """
  bcftools mpileup --threads ${task.cpus}\
  $known_sites_mpile_cmd \
  -Q ${params.min_base_quality}\
  -f $reffol/${refid}.fasta $bam | \
  bcftools call --threads ${task.cpus} \
  $call_varonly $known_sites_call_cmd \
  -O z -o ${outname}.vcf.gz
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
  // qual_filter = $params.min_snp_qual > 0 ? "  "
  """
  bcftools filter -e 'QUAL < $params.min_snp_qual' -O z -o ${name}.filtered.vcf.gz $vcf
  tabix ${name}.filtered.vcf.gz
  """
}

if (params.get_bed != false ){
  process getBed {
    tag "${name}"
    label 'env_medium'
    publishDir "${params.outdir}/bed", mode: "copy"

    input:
    set val(name), file(vcf), file(vcf_idx) from output_vcf

    output:
    set val(name), file("${name}.filtered.bed") into out_bed

    script:
    """
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%DP[\t%GT]\n" $vcf | awk '\$5 > $params.min_depth {print \$0}' > ${name}.filtered.bed
    """
  }
}


