#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process qc_post_alignment {
    clusterOptions = { 
        "--output=RNAseQC_${SRR}.out --error=RNAseQC_${SRR}.err" 
    }
    publishDir "${params.outDir}/QC_post_alignment", pattern: "*.metrics.tsv"
    publishDir "${params.outDir}/logs/RNAseQC", pattern: "*.{out,err}"
    
    input:
    tuple val(SRR), path(bam_file)

    output:
    path "*.metrics.tsv"
    path "*.{out,err}"

    script:
    """
    module load RNA-SeQC
    module load SAMtools

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: $SRR"

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Creating the index BAM file"

    samtools index $bam_file

    rnaseqc $params.gtfCollapsed $bam_file --sample $SRR --stranded rf --verbose .

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] RNAseQC of sample $SRR done"
    """
}