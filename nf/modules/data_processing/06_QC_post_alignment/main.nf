#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process collapsed_gtf {
    tag "Collapse_gtf_file"

    container "file://${projectDir}/containers/collapse_annotation.sif"

    input:
    path gtf

    output:
    path "annotation.collapsed_only.gtf"

    module "Apptainer"

    script:
    """
    input_gtf=$gtf
    filename=\$(basename $gtf)
    extension=\${filename##*.}
    if [[ \$extension == "gz" ]]; then
        gzip -d \$input_gtf > annotation.gtf
        input_gtf=annotation.gtf
    fi

    collapse_annotation.py \\
        --collapse_only \$input_gtf \\
        annotation.collapsed_only.gtf
    """
}

process qc_post_alignment {
    tag "$SRR"

    input:
    tuple val(SRR), path(bam_file)
    path collapse_gtf

    output:
    path "*.metrics.tsv", emit: "rnaseqc_metrics"

    module "RNA-SeQC:SAMtools"

    script:
    """
    samtools index $bam_file
    rnaseqc $collapse_gtf $bam_file --sample $SRR --stranded rf --verbose .
    """
}