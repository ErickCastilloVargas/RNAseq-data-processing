#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process fastQC_post_trimming {
    // Define inputs:
    input:
    path sample_dir from channel.fromPath("${params.out_dir}/trimmomatic_trimming/*/")

    // Define output: FastQC outputs post trimming
    output:
    path "${params.out_dir}/fastQC_reports_post_trimming/*"

    // Define the script
    script:
    """
    # Create the log dir 
    mkdir -p ${params.out_dir}/logs/fastQC_post_trimming
    
    # Load the module for fastQC
    module load FastQC

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: ${sample_dir.name}"

    # Output directory
    out_dir="${params.out_dir}/fastQC_reports_post_trimming"
    mkdir -p $out_dir

    # Input FASTQ files
    fastq1="${sample_dir}/${sample_dir.name}/${sample_dir.name}_1_paired.fastq.gz"
    fastq2="${sample_dir}/${sample_dir.name}/${sample_dir.name}_2_paired.fastq.gz"

    fastqc -t ${params.fastQC.threads} \
        ${fastq1} \
        ${fastq2} \
        -o ${out_dir}/

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] FastQC post trimming of sample ${sample_dir.name} done"
    """
}

process multiQC_post_trimming {
    // Define inputs:
    input:
    path fastQC_reports from channel.fromPath("${params.out_dir}/fastQC_reports_post_trimming")

    // Define output: FastQC outputs post trimming
    output:
    path "${params.out_dir}/multiQC_reports_post_trimming/*"

    // Define the script
    script:
    """
    # Create the log dir
    mkdir -p ${params.out_dir}/logs/multiQC_post_trimming

    # Load the required module
    module load MultiQC

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting multiQC post trimming ..."

    # Output dir for the multiQC reports
    out_dir="${params.out_dir}/multiQC_reports_post_trimming"
    mkdir -p "$out_dir"  

    # MultiQC for forward and reverse reads (they are paired-end)
    multiqc ${fastQC_reports}/*_1_paired_fastqc.zip -o ${out_dir} -n reports_1_forward_reads
    multiqc ${fastQC_reports}/*_2_paired_fastqc.zip -o ${out_dir} -n reports_2_reverse_reads

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] MultiQC post trimming done"
    """
}