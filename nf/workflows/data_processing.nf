#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include all modules
include { build_salmon_index; salmon_quant; infer_strandness; merge_strandness_results; mapping_rate; merge_mapping_rates } from "../modules/data_processing/00_Salmon"
include { adapters_trimming } from "../modules/data_processing/01_adapter_trimming"
include { fastQC_post_trimming } from "../modules/data_processing/02_QC_post_trimming"
include { star_and_RSEM_index_building } from "../modules/data_processing/03_STAR_and_RSEM_index_building"
include { star_alignment } from "../modules/data_processing/04_STAR_alignment"
include { samtools_sort_by_coordinates } from "../modules/data_processing/05_Sort_by_coordinates"
include { qc_post_alignment; collapsed_gtf } from "../modules/data_processing/06_QC_post_alignment"
include { rsem_transcript_quantification } from "../modules/data_processing/07_RSEM_transcript_quantification"
include { get_alignment_files } from "../modules/data_processing/get_alignment_files"
include { multiQC_final_report; agregate_rnaseqc_metrics } from "../modules/data_processing/final_QC_report"

workflow {
    def raw_samples_ch = Channel.fromFilePairs("${params.sampleDir}/SRR*_{1,2}.fastq.gz")
    
    def trimmed_samples_ch = Channel.empty()
    def fastqc_trimmed_ch = Channel.empty()
    if (params.adapter_trimming) {
        // Adapter Trimming
        def adapters_ch = Channel.fromPath(params.adapters)

        adapters_trimming(raw_samples_ch, adapters_ch.first())
        trimmed_samples_ch = adapters_trimming.out.trimmed_files

        // QC after trimming
        fastQC_post_trimming(trimmed_samples_ch)
        fastqc_trimmed_ch = fastQC_post_trimming.out.zip_files.collect()
    }

    // Salmon quantification (if pseudo aligner with salmon is specified) and strandness inference
    def transcriptome_ch = Channel.fromPath(params.transcriptome)
    def salmon_input_ch = params.adapter_trimming ? trimmed_samples_ch : raw_samples_ch

    build_salmon_index(transcriptome_ch)
    salmon_quant(build_salmon_index.out.first(), salmon_input_ch)
    infer_strandness(salmon_quant.out.lib_format_counts)
    strandness_SRR_ch = infer_strandness.out.map{ tuple -> tuple[0] }
    strandness_files_ch = infer_strandness.out.map{ tuple -> tuple[1] }
    merge_strandness_results(strandness_SRR_ch.collect(), strandness_files_ch.collect())

    mapping_rate_info_ch = Channel.empty()
    if (params.pseudoAligner == "Salmon") {
        mapping_rate(salmon_quant.out.mapping_rate_info)
        mapping_rate_SRR_ch = mapping_rate.out.map{ tuple -> tuple[0] }
        mapping_rate_files_ch = mapping_rate.out.map{ tuple -> tuple[1] }
        mapping_rate_info_ch = merge_mapping_rates(mapping_rate_SRR_ch.collect(), mapping_rate_files_ch.collect())
    }

    // STAR alignment and RSEM quantification if alignment is not skipped
    def rnaseq_metrics_ch = Channel.empty()
    if (!params.skipAlignment) {
        // STAR and RSEM index building
        def genome_ch = Channel.fromPath(params.genome)
        def gtf_ch = Channel.fromPath(params.gtf)

        star_and_RSEM_index_building(genome_ch, gtf_ch)

        // STAR Alignment
        def star_input_ch = params.adapter_trimming ? trimmed_samples_ch : raw_samples_ch
        star_alignment(star_input_ch, star_and_RSEM_index_building.out.star_index.first())

        // SAMtools sort by coordinates
        samtools_sort_by_coordinates(star_alignment.out.unsorted_bam_files)

        if (params.getAlignments) {
            get_alignment_files(samtools_sort_by_coordinates.out.sorted_bam_files, genome_ch.first())
        }

        // QC after alignment
        collapsed_gtf(gtf_ch)
        qc_post_alignment(samtools_sort_by_coordinates.out.sorted_bam_files, collapsed_gtf.out.first())

        // RSEM Transcript Quantification
        rsem_transcript_quantification(star_alignment.out.tr_bam_files, star_and_RSEM_index_building.out.rsem_index.first())

        // For the final QC report
        rnaseq_metrics_ch = agregate_rnaseqc_metrics(qc_post_alignment.out.rnaseqc_metrics.collect())
    }

    // Final QC report
    multiQC_final_report(rnaseq_metrics_ch, fastqc_trimmed_ch, merge_strandness_results.out, mapping_rate_info_ch)
}