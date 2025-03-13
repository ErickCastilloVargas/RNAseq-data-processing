#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
    // Include the modules 
    include { fastQC_raw_data; multiQC_raw_data} from "./modules/raw_data_QC"

    // Perform FastQc and multiQC of the raw data
    fastQC_raw_data()
    multiQC_raw_data()
}