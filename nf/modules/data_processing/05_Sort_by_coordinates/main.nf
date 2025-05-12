#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process samtools_sort_by_coordinates {
    clusterOptions = {
        "--output=SAMtools_${SRR}.out --error=SAMtools_${SRR}.err"
    }
    publishDir "${params.outDir}/SAMtools_sort_by_coordinates", pattern: "*sortByCoordinates.bam"
    publishDir "${params.outDir}/logs/SAMtools", pattern: "*.{out,err}"

    input:
    tuple val(SRR), path(bam_file)

    output:
    tuple val(SRR), path("*sortByCoordinates.bam"), emit: sorted_bam_files
    path("*.{out,err}")

    script:
    """
    module load SAMtools

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: $SRR"

    samtools sort -@ $params.samtools_threads -m 4G -o ${SRR}_sortByCoordinates.bam $bam_file

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Sample $SRR sorted by coordinates"
    """

}



