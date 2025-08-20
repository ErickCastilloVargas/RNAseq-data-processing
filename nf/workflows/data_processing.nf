#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include all modules
include { build_salmon_index; infer_strandness; merge_strandness_results } from "../modules/data_processing/00_infer_strandness"
include { adapters_trimming } from "../modules/data_processing/01_adapter_trimming"
include { fastQC_post_trimming } from "../modules/data_processing/02_QC_post_trimming"
include { star_and_RSEM_index_building } from "../modules/data_processing/03_STAR_and_RSEM_index_building"
include { star_alignment } from "../modules/data_processing/04_STAR_alignment"
include { samtools_sort_by_coordinates } from "../modules/data_processing/05_Sort_by_coordinates"
include { qc_post_alignment; collapsed_gtf } from "../modules/data_processing/06_QC_post_alignment"
include { rsem_transcript_quantification } from "../modules/data_processing/07_RSEM_transcript_quantification"
include { get_bam_files } from "../modules/data_processing/get_BAM_files"
include { multiQC_report_post_alignment; agregate_rnaseqc_metrics } from "../modules/data_processing/report_post_alignment"

workflow {
    def raw_samples_ch = Channel.fromFilePairs("${params.sampleDir}/SRR*_{1,2}.fastq.gz")

    // Infer strandness
    def transcriptome_ch = Channel.fromPath(params.transcriptome)

    build_salmon_index(transcriptome_ch)
    infer_strandness(build_salmon_index.out.first(), raw_samples_ch)
    strandness_SRR_ch = infer_strandness.out.map{ tuple -> tuple[0] }
    strandness_files_ch = infer_strandness.out.map{ tuple -> tuple[1] }
    merge_strandness_results(strandness_SRR_ch.collect(), strandness_files_ch.collect())

    // STAR and RSEM index building
    def genome_ch = Channel.fromPath(params.genome)
    def gtf_ch = Channel.fromPath(params.gtf)

    star_and_RSEM_index_building(genome_ch, gtf_ch)

    def trimmed_samples_ch = Channel.empty()
    if (params.adapter_trimming) {
        // Adapter Trimming
        def adapters_ch = Channel.fromPath(params.adapters)

        adapters_trimming(raw_samples_ch, adapters_ch.first())
        trimmed_samples_ch = adapters_trimming.out.trimmed_files

        // QC after trimming
        fastQC_post_trimming(trimmed_samples_ch)
    }

    // STAR Alignment
    def star_input_ch = params.adapter_trimming ? trimmed_samples_ch : raw_samples_ch
    star_alignment(star_input_ch, star_and_RSEM_index_building.out.star_index.first())

    // SAMtools sort by coordinates
    samtools_sort_by_coordinates(star_alignment.out.unsorted_bam_files)

    if (params.getBAM) {
        get_bam_files(samtools_sort_by_coordinates.out.sorted_bam_files, genome_ch.first())
    }

    // QC after alignment
    collapsed_gtf(gtf_ch)
    qc_post_alignment(samtools_sort_by_coordinates.out.sorted_bam_files, collapsed_gtf.out.first())

    // RSEM Transcript Quantification
    rsem_transcript_quantification(star_alignment.out.tr_bam_files, star_and_RSEM_index_building.out.rsem_index.first())

    // Final QC report
    agregate_rnaseqc_metrics(qc_post_alignment.out.rnaseqc_metrics.collect())
    multiQC_report_post_alignment(agregate_rnaseqc_metrics.out, fastQC_post_trimming.out.zip_files.collect(), merge_strandness_results.out)
}