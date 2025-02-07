nextflow.enable.dsl=2

process star_and_RSEM_index_building {
    // Define input: List of directories (sample directories) within the trimmomatic output directory
    input:
    path reference_genome from params.index_building.reference_genome
    path gtf_file from params.index_building.gtf_file
    val threads from params.index_building.threads

    // Define output: FastQC outputs post trimming
    output:
    path "${params.index_building.out_dir_STAR_index}/*"
    path "${params.index_building.out_dir_RSEM_index}/*"

    // Define the script
    script:
    """
    # STAR index generation
    STAR --runMode genomeGenerate \
        --genomeDir ${params.index_building.out_dir_STAR_index} \
        --genomeFastaFiles $reference_genome \
        --sjdbGTFfile $gtf_file \
        --sjdbOverhang 100 \
        --runThreadN $threads
     
    echo "STAR index created"

    # RSEM index generation
    rsem-prepare-reference --num-threads $threads \
        --gtf $gtf_file \
        $reference_genome \
        ${params.index_building.out_dir_RSEM_index}/rsem_reference
    
    echo "RSEM index created"
    """
}