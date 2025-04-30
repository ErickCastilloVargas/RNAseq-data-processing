#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process fastQC_post_trimming {
    clusterOptions = { 
        "--output=post_trimming_fastQC_${SRR}.out --error=post_trimming_fastQC_${SRR}.err"
    }
    publishDir "results/fastQC_reports_post_trimming", pattern: "*.{html,zip}"
    publishDir "results/logs/fastQC_post_trimming", pattern: "*.{out,err}"
    
    input:
    tuple val(SRR), path(fastq_files)

    output:
    path "*.zip", emit: "zip_files"
    path "*.html"
    path "*.{out,err}"

    script:
    """
    module load FastQC
    
    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: ${SRR}"
    
    fastqc -t ${params.fastQC.threads} \\
        $fastq_files \\
        -o .
    
    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] FastQC of post trimming sample ${SRR} done"
    """
}

process multiQC_post_trimming {
    clusterOptions = "--output=post_trimming_multiQC.out --error=post_trimming_multiQC.err"

    publishDir "results/multiQC_reports_post_trimming", pattern: "*.html"
    publishDir "results/logs/multiQC_post_trimming", pattern: "*.{out,err}"

    input:
    path fastQC_reports

    output:
    path "*.html"
    path "*.{out,err}"

    script:
    """
    module load MultiQC

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Starting multiQC of post_trimming data ..."

    multiqc $fastQC_reports -n multiQC_report_post_trimming

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] MultiQC of post_trimming data done"
    """
}