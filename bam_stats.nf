/*
Nextflow script for calculating statistics on the bam file
*/

params.input_bam = false
params.outdir = false


ch_input_bam = Channel
        .fromPath( params.input_bam )
        .map{ it -> [file(it).baseName, file(it)] }

custom_runName = params.run_name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

process bamStat {
    tag "$name"
    label 'env_medium'
    publishDir "${params.outdir}/${name}/", mode: 'copy'

    input:
    set val(name), file(bam) from ch_input_bam

    output:
    file "${name}_flagstat_report.txt" into ch_flagstat_results_for_multiqc
    file "${name}_stats_report.txt" into ch_samtools_stats_results_for_multiqc
    file "${name}_qualimap" into ch_qualimap_results_for_multiqc

    script:
    // samtools view -b -o ${name}.bam -S $sam
    // samtools sort -m 5G --threads ${task.cpus} -o ${name}.sorted.bam ${name}.bam
    """
    samtools index $bam
    samtools flagstat $bam > ${name}_flagstat_report.txt
    samtools stats ${bam} > ${name}_stats_report.txt

    qualimap bamqc -bam ${bam} -outdir ${name}_qualimap --collect-overlap-pairs --java-mem-size=10G -nt ${task.cpus}
    """
}

ch_config_for_multiqc = Channel
    .fromPath(params.multiqc_config, checkIfExists: true)
    .ifEmpty { exit 1, "multiqc config file not found: ${params.multiqc_config}" }


process multiQC {
    label 'env_medium'
    tag "${params.outdir}/MultiQC/$ofilename"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from ch_config_for_multiqc
    file ('qualimap/*') from ch_qualimap_results_for_multiqc.collect().ifEmpty([])
    file ('samtools/*') from ch_flagstat_results_for_multiqc.flatten().collect().ifEmpty([])
    file ('samtools/*') from ch_samtools_stats_results_for_multiqc.flatten().collect().ifEmpty([])

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