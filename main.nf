#!/usr/bin/env nextflow
/*
========================================================================================
             GATK HaplotypeCaller B E S T - P R A C T I C E
========================================================================================
 GATK 4 HaplotypeCaller Best Practice Analysis Pipeline.
 https://gatkforums.broadinstitute.org/gatk/discussion/11145/germline-short-variant-discovery-snps-indels
 #### Homepage / Documentation

 #### Authors
 Rahul Pisupati <rahul.pisupati@gmi.oeaw.ac.at>
----------------------------------------------------------------------------------------
*/


/*
Simply run this

nextflow run main.nf --reads "*bam" --file_ext bam --fasta ~/TAiR10_ARABIDOPSIS/TAIR10_wholeGenome.fasta --outdir output_folder

Make sure there is no hash tag in the file name. GATK 4 (GenomicsDBImport) isnt supporting this yet.
*/

/*
 * SET UP CONFIGURATION VARIABLES
 */
params.known_sites = false  // Known sites (VCF file) DB
// Store the chromosomes in a channel for easier workload scattering on large cohort
chromosomes = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]

params.project = "the1001genomes"
params.outdir = './snpcall'
params.cohort = "allsample"
params.run_name = false

params.fasta = false   // reference fasta file
params.file_ext = "bam"  // please change this accordingly..
params.singleEnd = false


params.clusterOptions = false
params.email = false
params.plaintext_email = false

params.saveTrimmed = false
params.notrim = false   // we trim the reads by default
params.illumina = true
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
  if( ! nextflow.version.matches(">= $params.nf_required_version") ){
    throw GroovyException('Nextflow version too old')
  }
} catch (all) {
  log.error "====================================================\n" +
            "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}

/*
 * Create a channel for input read files
 */
build_index = false
if ( params.fasta ){
  genome = file(params.fasta)
  reffol = genome.parent
  refid = genome.baseName
  if( !genome.exists() ) exit 1, "Reference fasta file not found: ${params.fasta}"
  bwa_indices = Channel
    .fromPath( "$reffol/${refid}.fasta.b*" )
    .ifEmpty { build_index = true }
    .subscribe onComplete: { checked_genome_index = true }
} else {
  exit 1, "Provide reference fasta file. Ex., --fasta file_path"
}

num_files = 1
if ( params.file_ext == 'fastq' ){
  num_files = params.singleEnd ? 1 : 2
}
read_files_processing = Channel
    .fromFilePairs( params.reads, size: num_files )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }

if ( params.known_sites != false ){
	known_sites_vcf = Channel
		.fromPath( params.known_sites )
		.map{ it -> [it, file("$it" + ".*")] }
	perform_bqsr = true
} else {
	perform_bqsr = false
}

/// ______________________________

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.run_name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

log.info "=================================================="
log.info " nf-haplocaller : SNP calling Best Practice v${params.version}"
log.info "=================================================="
def summary = [:]
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Reads']          = params.reads
summary['Data Type']      = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']         = params.fasta
if(params.notrim)       summary['Trimming Step'] = 'Skipped'
summary['Save Trimmed']   = params.saveTrimmed ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current user']   = "$USER"
summary['Current home']   = "$HOME"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
// summary['Container']      = process.container
// summary['Conda environment'] = process.conda
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "========================================="


/*
* 1. Create a channel for checking bwa index for genome ref
*/
if (build_index == true){
  process makeBWAindex {
      publishDir "${reffol}", mode: 'copy'
      label 'env_small'

      input:
      file genome

      output:
      file "${refid}.fasta.*" into bwa_index
      file "${refid}.dict" into fasta_dict

      script:
      """
      samtools faidx ${genome}
      bwa index $genome
      picard  CreateSequenceDictionary -R $genome -O ${refid}.dict
      """
  }
} else {
  bwa_index = Channel
    .fromPath( "$reffol/${refid}.fasta.*" )
  fasta_dict = Channel
    .fromPath( "$reffol/${refid}.dict" )
}

/*
* 2. Generate FASTQ from BAM file
*/
if (params.file_ext == "fastq"){
  read_files_processing.into { read_files_fastqc; read_files_trimming }
} else {

  process extractFastq {
    tag "$name"
    storeDir "${workflow.workDir}/rawreads"
    label 'env_small'

    input:
    set val(name), file(reads) from read_files_processing

    output:
    set val(name), file("${name}*fastq") into files_fastq

    script:
    if (params.singleEnd) {
      if (reads.getExtension() == "sra") {
        """
        fastq-dump $reads
        """
      } else if (reads.getExtension() == "bam") {
        """
        picard SamToFastq\
        I=$reads FASTQ=${name}.fastq VALIDATION_STRINGENCY=LENIENT
        """
      }
    } else {
      if (reads[0].getExtension() == "sra") {
        """
        fastq-dump --split-files $reads
        """
      } else if (reads.getExtension() == "bam") {
        """
        picard SamToFastq \
        I=$reads FASTQ=${name}_1.fastq SECOND_END_FASTQ=${name}_2.fastq \
        VALIDATION_STRINGENCY=LENIENT
        """
      }
    }
  }
  files_fastq.into{ read_files_fastqc; read_files_trimming}

}

/*
* 3. FastQC for the input files
*/
process fastqc {
    tag "$name"
    label 'env_small'
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file '*_fastqc.{zip,html}' into ch_fastqc_results_for_multiqc

    script:
    """
    fastqc -q $reads
    """
}

/*
* 4. Trimming the reads
*/
if(params.notrim){
    read_files_trimming.set{ trimmed_reads }
    ch_trimgalore_fastqc_reports_for_multiqc = Channel.create()
} else {
    process trimReads {
        tag "$name"
        label 'env_medium'
        publishDir "${params.outdir}/qc/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }

        input:
        set val(name), file(reads) from read_files_trimming

        output:
        set val(name), file('*fq.gz') into trimmed_reads
        file "*trimming_report.txt" into trimgalore_results
        file "*_fastqc.{zip,html}" into ch_trimgalore_fastqc_reports_for_multiqc

        script:
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        illumina = params.illumina ? "--illumina" : ''
        if (params.singleEnd) {
            """
            trim_galore --fastqc --gzip $illumina  $c_r1 $tpc_r1 $reads
            """
        } else {
            """
            trim_galore --paired --fastqc --gzip $illumina $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
            """
        }
    }
}


/*
* 5. Aligning the reads -- BWA-MEM
*/
process alignReads {
  tag "$name"
  label 'env_large'

  input:
  set val(name), file(reads) from trimmed_reads
  file genome
  file indices from bwa_index.collect()
  file fa_dict from fasta_dict.collect()

  output:
  set val(name), file("${name}.sam") into aligned_sam

  script:
  """
  bwa mem -t ${task.cpus} $genome $reads > ${name}.sam
  """
}

/*
* 6. Processing sam to bam and sorting the bam
*/
process processBam {
  tag "$name"
  label 'env_medium'
  publishDir "${params.outdir}/qc/alignstats_samtools/", mode: 'copy',
  saveAs: {filename ->
    if (filename.indexOf(".txt") > 0) "$filename"
    else null
  }


  input:
  set val(name), file(sam) from aligned_sam

  output:
  set val(name), file("${name}.sorted.bam") into sorted_bam
  file "${name}_flagstat_report.txt" into ch_flagstat_results_for_multiqc
  file "${name}_stats_report.txt" into ch_samtools_stats_results_for_multiqc

  script:
  """
  samtools view -b -o ${name}.bam -S $sam
  samtools sort -m 5G --threads ${task.cpus} -o ${name}.sorted.bam ${name}.bam

  samtools index ${name}.sorted.bam
  samtools flagstat ${name}.sorted.bam > ${name}_flagstat_report.txt
  samtools stats ${name}.sorted.bam > ${name}_stats_report.txt
  """
}

/*
* 7. Picard tool on bam file to remove duplicates
*/
process picardBam {
  tag "$name"
  label 'env_small'
  publishDir "${params.outdir}/", mode: 'copy',
  saveAs: {filename ->
    if (filename.indexOf(".bam") > 0) "alignedBam/$filename"
    else if (filename.indexOf(".bai") > 0) "alignedBam/$filename"
    else "qc/markdup_picard/$filename"
  }

  input:
  set val(name), file(bam) from sorted_bam

  output:
  set val(name), file("${name}.sorted.mkdup.bam"), file("${name}.sorted.mkdup.bam.bai") into (modified_bam, modified_bam_for_quali)
  file "${name}.metrics.txt" into ch_markDups_results_for_multiqc

  script:
  """
  picard MarkDuplicates\
  I=$bam O=${name}.dedup.bam METRICS_FILE=${name}.metrics.txt
  picard AddOrReplaceReadGroups\
  I=${name}.dedup.bam O=${name}.sorted.mkdup.bam\
  ID=$name LB=$name PL=illumina PU=none SM=$name
  samtools index ${name}.sorted.mkdup.bam
  """
}


/*
* 7.1 get qualimap for results on alignedBam
*/

process qcBam {
    tag "$name"
    label 'env_small'
    publishDir "${params.outdir}/qc/bamstats_qualimap", mode: 'copy'

    input:
    set val(name), file(bam), file(bam_index) from modified_bam_for_quali

    output:
    file "${name}_qualimap" into ch_qualimap_results_for_multiqc

    script:
    """
    qualimap bamqc \\
        -bam ${bam} \\
        -outdir ${name}_qualimap \\
        --collect-overlap-pairs \\
        --java-mem-size=${task.memory.toGiga()}G \\
        -nt ${task.cpus}
    """
}

/*
* 8. GATK to realign the reads at the positions where there are indels

This has been deprecated in GATK4

process realignBam {
  tag "$name"
  label 'env_gatk_small'

  input:
  set val(name), file(bam) from modified_bam
  set val(name), file(bam_index) from modified_bam_index

  output:
  set val(name), file("${name}.realignedBam.bam") into realigned_bam
  set val(name), file("${name}.realignedBam.bam.bai") into realigned_bam_index


  script:
  """
  gatk --java-options '-Djava.io.tmpdir=${params.tmpdir}' RealignerTargetCreator\
  -R $reffol/${refid}.fasta\
  -I $bam -o ${name}.intervals

  gatk --java-options '-Djava.io.tmpdir=${params.tmpdir}' IndelRealigner\
  -R $reffol/${refid}.fasta\
  -I $bam -targetIntervals ${name}.intervals\
  -o ${name}.realignedBam.bam

  samtools index ${name}.realignedBam.bam
  """
}
*/

/*
* 8.1 Base recalibration here
*/

if (perform_bqsr == true){
  process baseRecall {
    tag "$name"
    publishDir "${params.outdir}/", mode: 'copy',
    saveAs: {filename ->
      if (filename.indexOf("bam") > 0) "recallBam/$filename"
      else if (filename.indexOf("bai") > 0) "recallBam/$filename"
      else "qc/baseRecallstats_bqsr/$filename"
    }
    label 'env_large'

    input:
    set val(name), file(bam), file(bam_index) from modified_bam
    set file(vcf_db), file(vcf_db_idx) from known_sites_vcf.collect()

    output:
    set val(name), file("${name}.sorted.mkdup.recall.bam"), file("${name}.sorted.mkdup.recall.bam.bai") into recall_bam
    file("${bam}.BQSR.pdf") into stats_bqsr
    file ("${name}.recall_data.table") into ch_baserecal_results_for_multiqc

    script:
    """
    gatk BaseRecalibrator\
    --tmp-dir=${params.tmpdir} \
    -R $reffol/${refid}.fasta\
    -I $bam --known-sites $vcf_db\
    -O ${name}.recall_data.table

    gatk ApplyBQSR\
    --tmp-dir=${params.tmpdir} \
    -R $reffol/${refid}.fasta\
    -I $bam --bqsr-recal-file ${name}.recall_data.table\
    -O ${name}.sorted.mkdup.recall.bam

    gatk BaseRecalibrator\
    --tmp-dir=${params.tmpdir} \
    -R $reffol/${refid}.fasta\
    -I ${name}.sorted.mkdup.recall.bam --known-sites $vcf_db\
    -O ${name}.after_recal_data.table

    gatk AnalyzeCovariates\
    --tmp-dir=${params.tmpdir} \
    -before ${name}.recall_data.table -after ${name}.after_recal_data.table\
    -plots ${name}.BQSR.pdf

    samtools index ${name}.sorted.mkdup.recall.bam
    """
  }
} else {
  modified_bam.set{ recall_bam }
  ch_baserecal_results_for_multiqc = Channel.from([])
}


/*
* 9. GATK HaplotypeCaller for the SNPs
*/
process doSNPcall {
  tag "$name"
  publishDir "${params.outdir}/sampleGVCF", mode: 'copy'
  label 'env_large'

  input:
  set val(name), file(bam), file(bam_index) from recall_bam

  output:
  file("${name}.g.vcf.gz") into raw_gvcf
  file("${name}.g.vcf.gz.tbi") into raw_gvcf_index

  script:
  """
  gatk  HaplotypeCaller\
  --tmp-dir ${params.tmpdir}\
  -R $reffol/${refid}.fasta \
  -I $bam --output ${name}.g.vcf.gz \
  -ERC GVCF \
  """
  //-variant_index_type LINEAR -variant_index_parameter 128000
}


/*
* 10.1 GenomicsDBImport for all the files
*/
chromosomes_ch = Channel.from( chromosomes )

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
  set file("${name}.BIALLELIC.SNPs.vcf"), file("${name}.BIALLELIC.SNPs.vcf.idx") into filter_vcf

  script:
  """
  gatk SelectVariants --java-options "-Xmx24g -Xms24g" \
    --tmp-dir ${params.tmpdir}\
    -R $reffol/${refid}.fasta \
    -V $gvcf \
    -select-type SNP --restrict-alleles-to BIALLELIC \
    -se "${name}" \
    -O ${name}.BIALLELIC.SNPs.vcf
  """
}

/*
 * STEP 000 - MultiQC
 */
ch_config_for_multiqc = Channel
    .fromPath(params.multiqc_config, checkIfExists: true)
    .ifEmpty { exit 1, "multiqc config file not found: ${params.multiqc_config}" }

/*
 * Parse software version numbers
 */
// process get_software_versions {
//   tag 'env_small'
//
//   output:
//   file 'software_versions_mqc.yaml' into ch_software_versions_yaml_for_multiqc
//
//   script:
//   """
//   echo "$workflow.manifest.version" &> v_ngi_methylseq.txt
//   echo "$workflow.nextflow.version" &> v_nextflow.txt
//   fastqc --version &> v_fastqc.txt
//   trim_galore --version &> v_trim_galore.txt
//   bwa &> v_bwa.txt 2>&1 || true
//   picard MarkDuplicates --version &> v_picard_markdups.txt 2>&1 || true
//   gatk --version &> v_gatk.txt
//   multiqc --version &> v_multiqc.txt
//   """
// }

process multiQC {
    label 'env_medium'
    tag "${params.outdir}/MultiQC/$ofilename"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from ch_config_for_multiqc
    file ('fastqc/*') from ch_fastqc_results_for_multiqc.collect().ifEmpty([])
    file ('trimgalore/*') from ch_trimgalore_fastqc_reports_for_multiqc.collect().ifEmpty([])
    file ('picard/*') from ch_markDups_results_for_multiqc.flatten().collect().ifEmpty([])
    file ('qualimap/*') from ch_qualimap_results_for_multiqc.collect().ifEmpty([])
    file ('gatk/*') from ch_baserecal_results_for_multiqc.collect().ifEmpty([])
    file ('samtools/*') from ch_flagstat_results_for_multiqc.flatten().collect().ifEmpty([])
    file ('samtools/*') from ch_samtools_stats_results_for_multiqc.flatten().collect().ifEmpty([])
    // file ('software_versions/*') from ch_software_versions_yaml_for_multiqc.collect().ifEmpty([])

    output:
    file "*_report.html" into ch_multiqc_report
    file "*_data"
    file '.command.err'

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    if(custom_runName){
      rfilename = "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
      ofilename = rfilename+'.html'
    } else {
      rfilename = ''
      ofilename = 'multiqc_report.html'
    }
    """
    multiqc --ignore *_R2* -f $rtitle $rfilename --config $multiqc_config .
    """
}
