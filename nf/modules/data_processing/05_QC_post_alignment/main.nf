#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process qc_post_alignment {
    clusterOptions = { 
        "--output=RNAseQC_${SRR}.out --error=RNAseQC_${SRR}.err" 
    }
    publishDir "results/QC_post_alignment", pattern: "*.metrics.tsv"
    publishDir "results/logs/RNAseQC", pattern: "*.{out,err}"
    
    input:
    tuple val(SRR), path(bam_file)

    output:
    path "*.metrics.tsv"
    path "*.{out,err}"

    script:
    """
    module load RNA-SeQC

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: ${SRR}"

    rnaseqc ${params.rnaseQC.collapsed_gtf_file} ${bam_file} --sample ${SRR} --stranded rf --verbose .

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] RNAseQC of sample ${SRR} done"
    """
}