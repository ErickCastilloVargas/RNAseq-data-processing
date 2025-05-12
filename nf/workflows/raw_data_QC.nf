#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include the modules 
include { fastQC_raw_data; multiQC_raw_data} from "../modules/raw_data_QC/main.nf"

workflow {
    // Perform FastQC and multiQC of the raw data
    sample_ch = channel.fromFilePairs("${params.sampleDir}/SRR*_{1,2}.fastq.gz")
    fastQC_raw_data(sample_ch)

    multiQC_raw_data(fastQC_raw_data.out.zip_files.collect())
}