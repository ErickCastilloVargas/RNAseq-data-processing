#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process rsem_transcript_quantification {
    tag "$SRR"

    publishDir "${params.outDir}/Gene_level_quant", mode:"move", pattern: "*.genes.results"

    input:
    tuple val(SRR), path(bam_file)
    path rsem_index

    output:
    path "*.genes.results"

    module "RSEM"

    script:
    """
    rsem-calculate-expression \\
        --num-threads $params.rsem_threads \\
        --fragment-length-max 1000 \\
        --append-names \\
        --no-bam-output \\
        --paired-end \\
        --estimate-rspd \\
        --alignments \\
        $bam_file \\
        ${rsem_index}/hg38 \\
        $SRR
    """
}