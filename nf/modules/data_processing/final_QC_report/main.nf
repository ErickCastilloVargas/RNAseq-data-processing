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

    publishDir "${params.outDir}/Reports", mode: "move"

    input:
    path rnaseqc_metric_file
    path fastqc_reports
    path strandness
    path mapping_rate_info

    output:
    path "QC_report_final.html"

    module "MultiQC"

    script:
    """
    cp "${projectDir}/resources/multiqc_config.yaml" custom_config.yaml

    # Dynamic filtering of the config.yaml
    if [[ -z "$rnaseqc_metric_file" ]]; then
        # Remove section RNAseQC_metrics
        awk '/^RNAseQC_metrics:/ {skip=1; next} skip && /^[^ ]/ {skip=0} !skip' custom_config.yaml > tmp && mv tmp custom_config.yaml
    fi

    if [[ -z "$mapping_rate_info" ]]; then
        awk '/^mapping_rate_info:/ {skip=1; next} skip && /^[^ ]/ {skip=0} !skip' custom_config.yaml > tmp && mv tmp custom_config.yaml
    fi

    Arg=()
    [[ -n "$rnaseqc_metric_file" ]] && Arg+=("$rnaseqc_metric_file")
    [[ -n "$fastqc_reports" ]] && Arg+=("$fastqc_reports")
    [[ -n "$strandness" ]] && Arg+=("$strandness")
    [[ -n "$mapping_rate_info" ]] && Arg+=("$mapping_rate_info")

    multiqc "\${Arg[@]}" \
        -c custom_config.yaml \
        -n QC_report_final \
        -o .
    """
}
