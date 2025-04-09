#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process adapters_trimming {
    clusterOptions = { 
        "--cpus-per-task=${params.trimmomatic.threads} --output=trimmomatic_trimming_${SRR}.out --error=trimmomatic_trimming_${SRR}.err" 
    }
    publishDir "results/trimmomatic_adapter_trimming", pattern: "*_{1,2}_paired.fastq.gz"
    publishDir "results/logs/adapter_trimming", pattern: "*.{out,err}"
    
    input:
    tuple val(SRR), path(fastq_files)

    output:
    tuple val(SRR), path("*_{1,2}_paired.fastq.gz"), emit: "trimmed_files"
    path "*.{out,err}"

    script:
    """
    module load Trimmomatic

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: ${SRR}"

    java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads ${params.trimmomatic.threads} \\
        $fastq_files \\
        ${SRR}_1_paired.fastq.gz \\
        ${SRR}_1_unpaired.fastq.gz \\
        ${SRR}_2_paired.fastq.gz \\
        ${SRR}_2_unpaired.fastq.gz \\
        ILLUMINACLIP:${params.trimmomatic.adapters_file}:2:30:10 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] All adapters are removed from ${SRR}"
    """
}