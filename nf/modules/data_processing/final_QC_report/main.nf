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

process multiQC_final_report {
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
    Arg=""
    if [[ -n "$rnaseqc_metric_file" && -n "$fastqc_reports" ]]; then
        Arg="${rnaseqc_metric_file} ${fastqc_reports} ${strandness}"
        cp ${projectDir}/resources/multiqc_config.yaml custom_config.yaml
    elif [[ -n "$rnaseqc_metric_file" && -z "$fastqc_reports" ]]; then
        Arg="${rnaseqc_metric_file} ${strandness}"
        cp ${projectDir}/resources/multiqc_config.yaml custom_config.yaml
    elif [[ -z "$rnaseqc_metric_file" && -n "$fastqc_reports" ]]; then
        Arg="${fastqc_reports} ${strandness}"
        awk '/^RNAseQC_metrics:/ {skip=1; next} 
            skip && /^[^ ]/ {skip=0} 
            !skip' ${projectDir}/resources/multiqc_config.yaml > custom_config.yaml
    else 
        Arg="${strandness}"
        awk '/^RNAseQC_metrics:/ {skip=1; next} 
            skip && /^[^ ]/ {skip=0} 
            !skip' ${projectDir}/resources/multiqc_config.yaml > custom_config.yaml
    fi

    multiqc \\
        \$Arg \\
        -c custom_config.yaml \\
        -n QC_report_final \\
        -o .
    """
}
