#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process fastQC_raw_data {
    tag "$SRR"

    input:
    tuple val(SRR), path(fastq_files)

    output:
    path "*.zip", emit: "zip_files"

    module "FastQC"

    script:
    """
    fastqc -t $params.fastQC_threads \\
        $fastq_files \\
        -o .
    
    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] FastQC of raw sample $SRR done"
    """
}

process multiQC_raw_data {
    publishDir "${params.outDir}/Reports/raw_data", mode:"copy", pattern: "*.html"

    input:
    path fastQC_reports

    output:
    path "*.html"

    module "MultiQC"

    script:
    """
    multiqc $fastQC_reports -n multiQC_report_raw_data
    """
}