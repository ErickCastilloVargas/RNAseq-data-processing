#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process adapters_trimming {
    tag "$SRR" 

    input:
    tuple val(SRR), path(fastq_files)
    path adapters

    output:
    tuple val(SRR), path("*_{1,2}_paired.fastq.gz"), emit: "trimmed_files"

    module "Trimmomatic"

    script:
    """
    java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $params.trimmomatic_threads \\
        $fastq_files \\
        ${SRR}_1_paired.fastq.gz \\
        ${SRR}_1_unpaired.fastq.gz \\
        ${SRR}_2_paired.fastq.gz \\
        ${SRR}_2_unpaired.fastq.gz \\
        ILLUMINACLIP:${adapters}:2:30:10 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}