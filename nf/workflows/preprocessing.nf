#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
    // Include the modules 
    include { adapters_trimming } from "../modules/preprocessing/01_adapter_trimming/adapter_trimming.nf"
    include { qc_post_trimming } from "../modules/preprocessing/02_QC_post_trimming/FastQC_post_trimming.nf"
    include { star_alignment } from "../modules/preprocessing/04_STAR_alignment/STAR_alignment.nf"
    include { star_alignment_no_adapter_trimming } from "../modules/preprocessing/04_STAR_alignment/STAR_alignment.nf"
    include { qc_post_alignment } from "../modules/preprocessing/05_QC_post_alignment/QC_post_alignment.nf"
    include { rsem_transcript_quantification } from "../modules/preprocessing/06_RSEM_transcript_quantification/RSEM_transcript_quantification.nf"

    if (params.adapter_trimming) {
        // Step 1: Adapter Trimming
        adapters_trimming()

        // Step 2: QC after trimming
        qc_post_trimming()

        // Step 3: STAR Alignment
        star_alignment()
    } else {
        // Step 3: STAR Alignment without adapter trimming
        star_alignment_no_adapter_trimming()
    }

    // Step 4: QC after alignment
    qc_post_alignment()

    // Step 5: RSEM Transcript Quantification
    rsem_transcript_quantification()
}