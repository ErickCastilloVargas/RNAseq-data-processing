#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include all modules
include { adapters_trimming } from "../modules/data_processing/01_adapter_trimming"
include { fastQC_post_trimming; multiQC_post_trimming} from "../modules/data_processing/02_QC_post_trimming"
include { star_and_RSEM_index_building } from "../modules/data_processing/03_STAR_and_RSEM_index_building"
include { star_alignment } from "../modules/data_processing/04_STAR_alignment"
include { samtools_sort_by_coordinates } from "../modules/data_processing/05_Sort_by_coordinates"
include { qc_post_alignment } from "../modules/data_processing/06_QC_post_alignment"
include { rsem_transcript_quantification } from "../modules/data_processing/07_RSEM_transcript_quantification"

workflow {
    // Step 1: STAR and RSEM index building
    star_and_RSEM_index_building()

    def raw_samples_ch = Channel.fromFilePairs("${params.sampleDir}/SRR*_{1,2}.fastq.gz")

    def trimmed_samples_ch = Channel.empty()
    if (params.adapter_trimming) {
        // Adapter Trimming
        adapters_trimming(raw_samples_ch)
        trimmed_samples_ch = adapters_trimming.out.trimmed_files

        // QC after trimming
        fastQC_post_trimming(trimmed_samples_ch)
        multiQC_post_trimming(fastQC_post_trimming.out.zip_files.collect())
    }

    // Step 2: STAR Alignment
    def star_input_ch = params.adapter_trimming ? trimmed_samples_ch : raw_samples_ch
    star_alignment(star_input_ch, star_and_RSEM_index_building.out.star_index)

    // Step 3: SAMtools sort by coordinates
    samtools_sort_by_coordinates(star_alignment.out.unsorted_bam_files)

    // Step 4: QC after alignment
    qc_post_alignment(samtools_sort_by_coordinates.out.sorted_bam_files)

    // Step 5: RSEM Transcript Quantification
    rsem_transcript_quantification(star_alignment.out.tr_bam_files, star_and_RSEM_index_building.out.rsem_index)
}