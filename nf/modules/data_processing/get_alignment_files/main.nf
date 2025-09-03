#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process get_alignment_files {
    tag "$SRR"

    publishDir "${params.outDir}/Alignment_files", mode: "move", pattern: "*.cram"

    when:
    params.getAlignments
    
    input:
    tuple val(SRR), path(bam_file)
    path genome

    output:
    path "${SRR}.cram"

    module "SAMtools"

    script:
    """
    samtools faidx $genome
    samtools view -T $genome -C -o ${SRR}.cram $bam_file
    """
}
