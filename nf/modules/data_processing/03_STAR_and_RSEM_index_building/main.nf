#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process star_and_RSEM_index_building {
    clusterOptions = "--output=index_building.out --error=index_building.err"

    publishDir "${params.outDir}/indexes", pattern: "STAR_index"
    publishDir "${params.outDir}/indexes", pattern: "RSEM_index"
    publishDir "${params.outDir}/logs/STAR_RSEM_index_building", pattern: "*.{out,err}"

    output:
    path "STAR_index", emit: "star_index"
    path "RSEM_index", emit: "rsem_index"
    path "*.{out,err}"

    script:
    """
    module load STAR
    module load RSEM

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Starting to create STAR index"

    STAR --runMode genomeGenerate \\
        --genomeDir STAR_index \\
        --genomeFastaFiles $params.genome \\
        --sjdbGTFfile $params.gtf \\
        --sjdbOverhang 100 \\
        --genomeChrBinNbits 18 \\
        --runThreadN $params.index_building_threads
     
    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] STAR index created"

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Starting to create RSEM index"
    
    mkdir -p RSEM_index
    rsem-prepare-reference --num-threads $params.index_building_threads \\
        --gtf $params.gtf \\
        $params.genome \\
        RSEM_index/hg38
    
    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] RSEM index created"
    """
}