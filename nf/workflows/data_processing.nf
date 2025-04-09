#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include all modules
// include { adapters_trimming } from "../modules/data_processing/01_adapter_trimming"
// include { fastQC_post_trimming; multiQC_post_trimming} from "../modules/data_processing/02_QC_post_trimming"
include { star_and_RSEM_index_building } from "../modules/data_processing/03_STAR_and_RSEM_index_building"
// include { star_alignment; star_alignment_no_adapter_trimming } from "../modules/data_processing/04_STAR_alignment"
// include { qc_post_alignment } from "../modules/data_processing/05_QC_post_alignment"
// include { rsem_transcript_quantification } from "../modules/data_processing/06_RSEM_transcript_quantification"

workflow {
    // Step 1: STAR and RSEM index building
    star_and_RSEM_index_building()

    if (params.adapter_trimming) {
        // Step 1: Adapter Trimming
        sample_ch = channel.fromFilePairs("${params.main_sample_dir}/SRR*_{1,2}.fastq.gz")
        sample_ch.view { "Samples tuple encontrado: ${it}" }
        adapters_trimming(sample_ch)

        // Step 2: QC after trimming
        fastQC_post_trimming(adapters_trimming.out.trimmed_files)
        multiQC_post_trimming(fastQC_post_trimming.out.zip_files.collect())

        // Step 4: STAR Alignment
        star_alignment(adapters_trimming.out.trimmed_files, star_and_RSEM_index_building.out.star_index)
    } else {
        // Step 4: STAR Alignment without adapter trimming
        sample_ch = channel.fromFilePairs("${params.main_sample_dir}/SRR*_{1,2}.fastq.gz")
        star_alignment_no_adapter_trimming(sample_ch, star_and_RSEM_index_building.out.star_index)
    }

    // // Step 5: QC after alignment
    // qc_post_alignment()

    // // Step 6: RSEM Transcript Quantification
    // rsem_transcript_quantification()
}