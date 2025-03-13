#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
    // Include the modules 
    include { adapters_trimming } from "./modules/data_processing/01_adapter_trimming"
    include { fastQC_post_trimming; multiQC_post_trimming} from "./modules/data_processing/02_QC_post_trimming"
    include { star_and_RSEM_index_building } from "./modules/data_processing/03_STAR_and_RSEM_index_building"
    include { star_alignment; star_alignment_no_adapter_trimming } from "./modules/data_processing/04_STAR_alignment"
    include { qc_post_alignment } from "./modules/data_processing/05_QC_post_alignment"
    include { rsem_transcript_quantification } from "./modules/data_processing/06_RSEM_transcript_quantification"

    // Step 1: STAR and RSEM index building
    star_and_RSEM_index_building()

    if (params.adapter_trimming) {
        // Step 1: Adapter Trimming
        adapters_trimming()

        // Step 2: QC after trimming
        fastQC_post_trimming()
        multiQC_post_trimming()

        // Step 4: STAR Alignment
        star_alignment()
    } else {
        // Step 4: STAR Alignment without adapter trimming
        star_alignment_no_adapter_trimming()
    }

    // Step 5: QC after alignment
    qc_post_alignment()

    // Step 6: RSEM Transcript Quantification
    rsem_transcript_quantification()
}