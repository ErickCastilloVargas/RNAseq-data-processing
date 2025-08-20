#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process agregate_rnaseqc_metrics {
    tag "agregate_RNAseQC_metrics_files"

    input:
    path rnaseqc_metrics

    output:
    path "rnaseqc_metrics_full.tsv"

    script:
    """
    cut -f1 ${rnaseqc_metrics[0]} | paste -s > rnaseqc_metrics_full.tsv
    
    for file in $rnaseqc_metrics; do
        cut -f2 \$file | paste -s >> rnaseqc_metrics_full.tsv
    done
    """
}

process multiQC_report_post_alignment {
    tag "MultiQC_report"

    publishDir "${params.outDir}/Reports" 

    input:
    path rnaseqc_metric_file
    path fastqc_reports
    path strandness

    output:
    path "QC_report_final.html"

    module "MultiQC"

    script:
    """
    multiqc \\
        $rnaseqc_metric_file \\
        $fastqc_reports \\
        $strandness \\
        -c ${projectDir}/resources/multiqc_config.yaml \\
        -n QC_report_final \\
        -o .
    """
}
