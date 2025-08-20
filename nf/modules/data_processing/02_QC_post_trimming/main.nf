#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process fastQC_post_trimming {
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
    """
}