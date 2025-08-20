#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process samtools_sort_by_coordinates {
    tag "$SRR"

    input:
    tuple val(SRR), path(bam_file)

    output:
    tuple val(SRR), path("*sortByCoordinates.bam"), emit: sorted_bam_files

    module "SAMtools"

    script:
    """
    samtools sort -@ $params.samtools_threads -m 4G -o ${SRR}_sortByCoordinates.bam $bam_file
    """

}



