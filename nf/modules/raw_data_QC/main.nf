#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process fastQC_raw_data {
    clusterOptions = { 
        "--cpus-per-task=${params.fastQC.threads} --output=raw_fastQC_${SRR}.out --error=raw_fastQC_${SRR}.err" 
    }

    publishDir "results/fastQC_reports_raw_data", pattern: "*.{html,zip}"
    publishDir "results/logs/fastQC_raw_data", pattern: "*.{out,err}"

    input:
    tuple val(SRR), path(fastq_files)

    output:
    //path "*.{zip,html}", emit: "zip_files"
    path "*.zip", emit: "zip_files"
    path "*.html"
    path "*.{out,err}"

    script:
    """
    module load FastQC
    
    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: ${SRR}"
    
    fastqc -t ${params.fastQC.threads} \\
        $fastq_files \\
        -o .
    
    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] FastQC of raw sample ${SRR} done"
    """
}

process multiQC_raw_data {
    clusterOptions = "--cpus-per-task=${params.fastQC.threads} --output=raw_data_multiQC.out --error=raw_data_multiQC.err"

    publishDir "results/multiQC_reports_raw_data", pattern: "*.html"
    publishDir "results/logs/multiQC_raw_data", pattern: "*.{out,err}"

    input:
    path fastQC_reports

    output:
    path "*.html"
    path "*.{out,err}"

    script:
    """
    module load MultiQC

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Starting multiQC of raw data ..."

    multiqc $fastQC_reports -n multiQC_report_raw_data

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] MultiQC of raw data done"
    """
}