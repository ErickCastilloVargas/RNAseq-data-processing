#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process star_and_RSEM_index_building {
    tag "Indexes_building"

    input:
    path genome
    path gtf
    
    output:
    path "STAR_index", emit: "star_index"
    path "RSEM_index", emit: "rsem_index"

    module "STAR:RSEM"
    script:
    """
    STAR --runMode genomeGenerate \\
        --genomeDir STAR_index \\
        --genomeFastaFiles $genome \\
        --sjdbGTFfile $gtf \\
        --sjdbOverhang 100 \\
        --genomeChrBinNbits 18 \\
        --runThreadN $params.index_building_threads
    
    mkdir -p RSEM_index
    rsem-prepare-reference --num-threads $params.index_building_threads \\
        --gtf $gtf \\
        $genome \\
        RSEM_index/hg38
    """
}