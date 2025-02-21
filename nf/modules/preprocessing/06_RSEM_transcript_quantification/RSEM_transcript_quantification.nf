#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process rsem_transcript_quantification {
    // Define inputs:
    input:
    path star_sample_dir from file("${params.out_dir}/STAR_alignment/*/")
    path rsem_index from ${params.out_dir}/rsem_index
    val threads from params.rsem.threads

    // Define output: RSEM transcript quantification
    output:
    path "${params.out_dir}/rsem_quantification/${star_sample_dir.name}/*"

    // Define the script
    script:
    """
    # Create the log dir
    mkdir -p ${params.out_dir}/logs/rsem_transcript_quantification

    # Load the module for RSEM
    module load RSEM

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: ${star_sample_dir.name}"

    # STAR .bam file aligned to transcriptome
    star_bam_file="${star_sample_dir}/${star_sample_dir.name}/${star_sample_dir.name}_Aligned.toTranscriptome.out.bam"

    # Create the output dir
    out_dir="${params.out_dir}/rsem_quantification/${star_sample_dir.name}"
    mkdir -p $out_dir

    rsem-calculate-expression \
        --num-threads ${threads} \
        --fragment-length-max 1000 \
        --append-names \
        --no-bam-output \
        --paired-end \
        --estimate-rspd \
        --alignments \
        ${star_bam_file} \
        ${rsem_index}/rsem_reference ${out_dir}/${star_sample_dir.name}

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Transcript quantification for sample ${star_sample_dir.name} finished"
    """
}