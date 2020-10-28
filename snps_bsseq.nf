/*
 to call SNPs using BisSNP in known sites 

nextflow run snpsnf --input "*bam" --known_sites 1001g_vcf --fasta ref_seq/TAIR10_wholeGenome.fasta  --outdir filter_vcfs
*/

params.project = "the1001genomes"
params.input = false
params.outdir = './snpcall'
params.fasta = false
params.cohort = false // file name for output combined vcf

/* Below are additional filtering steps, only work on cohort analysis
1. filter sites based on minor allele frequency (maf)
2. missing fraction of missing calls
*/
params.filter_maf = false
params.min_maf = 0.2  //
params.max_na = 0.3 //

params.snpcaller = "bcftools"
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


//  SNP calling starts now
input_bams = Channel
        .fromPath( "${params.input}" )
        .map{ it -> [ file(it).baseName, file(it), file("${it}.bai") ] }



// process MethylExtract {
//   tag "${name}"
//   label "env_medium"
//   publishDir "${params.outdir}/raw_variants", mode: "copy"

//   input:
//   set val(name), file(bam) from modify_bam

//   output:
//   file("${name}.filtered.snpvcf.bed") into bsseq_vcfbed

//   script:
//   outname = bam.baseName 
//   // flagW=99,147 flagC=83,163 // set for paired end
//   """
//   mkdir out_dir
//   ln -s -r $bam out_dir/$bam

//   MethylExtract.pl seq=~/TAiR10_ARABIDOPSIS/TAIR10_wholeGenome.fasta \
//   inDir=./ outDir=out_dir \
//   flagW=99,147 flagC=83,163 \  
//   varFraction=0.05 maxPval=0.01 \
//   LastIgnor=3 minDepthSNV=2 \

  
//   bcftools mpileup --threads ${task.cpus}\
//   $known_sites_mpile_cmd \
//   -Q ${params.min_base_quality}\
//   -f $reffol/${refid}.fasta $bam | \
//   bcftools call --threads ${task.cpus} \
//   $call_varonly $known_sites_call_cmd \
//   -O z -o ${outname}.vcf.gz
  
//   bcftools filter -T $known_sites_targets ${name}.vcf.gz | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%AD][\t%GT]\n" |\
//   awk 'length(\$3) == 1 && length(\$4) == 1 {print \$0}'  > ${name}.snpvcf.bed
    
//   python $workflow.projectDir/scripts/01.filter_bsseq_variants.py \
//   -i ${name}.snpvcf.bed -o ${name}.filtered.snpvcf.bed
//   """
// }


/*
Using BCFtools to call genotypes at known sites
*/

if (params.snpcaller == "bcftools"){

  process BCFcall {
    tag "${name}"
    label "env_medium"
    publishDir "${params.outdir}/variants_bcftools", mode: "copy"

    input:
    set val(name), file(bam), file(bam_index) from input_bams

    output:
    file("${name}.snp.vcf.gz*") into bsseq_snps
    file("${name}.qual_filtered.vcf.gz*") into bsseq_qual_snps 

    script:
    outname = bam.baseName 
    // flagW=99,147 flagC=83,163 // set for paired end
    known_sites_mpile_cmd = params.known_sites != false ? "-T ${params.known_sites}" : ''
    known_sites_call_cmd = params.known_sites != false ? "-C alleles -m -T ${params.known_sites}" : ''
    """
    bcftools mpileup --threads ${task.cpus}\
    $known_sites_mpile_cmd \
    -Q ${params.min_base_quality}\
    -f $reffol/${refid}.fasta $bam | \
    bcftools call --threads ${task.cpus} \
    $known_sites_call_cmd \
    -O z -o ${name}.raw.vcf.gz

    bcftools filter -e 'QUAL < $params.min_snp_qual' \
    -Oz ${name}.raw.vcf.gz | \
    bcftools view -V indels > ${name}.snps.vcf
    bgzip ${name}.snps.vcf && tabix ${name}.snps.vcf.gz
    
    bcftools filter -e ' REF == "G" && ALT == "A" ' ${name}.snps.vcf |\
    bcftools filter -e ' REF == "C" && ALT == "T" ' > ${name}.qual_filtered.vcf
    bgzip ${name}.qual_filtered.vcf && tabix ${name}.qual_filtered.vcf.gz
  
    """
  }

}

/*
Using GATK haplotyper caller below to call SNPs
Caution: Havent tested SNPs called with GATK yet.. ended with some errors

*/

if (params.snpcaller == "gatk"){
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
    // Using GATK
    tag "${name}"
    label "env_medium"
    publishDir "${params.outdir}/raw_variants", mode: "copy"

    input:
    set val(name), file(bam) from modify_bam

    output:
    file("${name}.qual_filtered.vcf.gz*") into bsseq_vcfbed

    script:
    known_sites_cmd = params.known_sites != false ? "--alleles ${params.known_sites} --force-call-filtered-alleles" : ''
    known_sites_bcf = params.known_sites != false ? "-T ${params.known_sites}" : ''
    // --alleles ${known_sites_vcf} --force-call-filtered-alleles \
    """
    samtools index $bam 

    gatk --java-options '-Djava.io.tmpdir=${params.tmpdir}' HaplotypeCaller\
    -R $reffol/${refid}.fasta\
    -I $bam \
    -O ${name}.raw.vcf.gz \
    $known_sites_cmd \
    --output-mode EMIT_ALL_CONFIDENT_SITES

    bcftools filter -e 'QUAL < $params.min_snp_qual' \
    -Oz ${name}.raw.vcf.gz | \
    bcftools view -V indels | \
    bcftools filter -e ' REF == "G" && ALT == "A" ' |\
    bcftools filter -e ' REF == "C" && ALT == "T" ' > ${name}.qual_filtered.vcf

    bgzip ${name}.qual_filtered.vcf && tabix ${name}.qual_filtered.vcf.gz
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
    label "env_medium"
    publishDir "${params.outdir}/", mode: "copy"

    input:
    file(vcf) from bsseq_qual_snps.collect()

    output:
    file("${params.cohort}.qual_filtered.vcf.gz") into mergedVCF

    script:
    """
    bcftools merge *.vcf.gz -Oz -o ${params.cohort}.qual_filtered.vcf.gz
    """
  }
  if (params.filter_maf) {
    process filterMAF {
      lable "env_medium"
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