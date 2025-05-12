#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process rsem_transcript_quantification {
    clusterOptions = { 
        "--output=RSEM_${SRR}.out --error=RSEM_${SRR}.err" 
    }
    publishDir "${params.outDir}/RSEM_transcript_quant", pattern: "*.genes.results"
    publishDir "${params.outDir}/logs/RSEM", pattern: "*.{out,err}"

    input:
    tuple val(SRR), path(bam_file)
    path rsem_index

    output:
    path "*.genes.results"
    path "*.{out,err}"

    // Define the script
    script:
    """
    module load RSEM

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: $SRR"

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

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Transcript quantification for sample $SRR finished"
    """
}