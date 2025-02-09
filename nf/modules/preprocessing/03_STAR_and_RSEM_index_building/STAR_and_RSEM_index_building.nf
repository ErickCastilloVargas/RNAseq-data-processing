nextflow.enable.dsl=2

process star_and_RSEM_index_building {
    // Define inputs:
    input:
    path reference_genome from params.index_building.reference_genome
    path gtf_file from params.index_building.gtf_file
    val threads from params.index_building.threads

    // Define outputs: indexes of STAR and RSEM
    output:
    path "${params.out_dir}/star_index/*"
    path "${params.out_dir}/rsem_index/*"

    // Define the script
    script:
    """
    # Load the module for STAR and RSEM
    module load STAR
    module load RSEM

    # STAR index generation
    STAR --runMode genomeGenerate \
        --genomeDir ${params.out_dir}/star_index \
        --genomeFastaFiles $reference_genome \
        --sjdbGTFfile $gtf_file \
        --sjdbOverhang 100 \
        --runThreadN $threads
     
    echo "STAR index created"

    # RSEM index generation
    rsem-prepare-reference --num-threads $threads \
        --gtf $gtf_file \
        $reference_genome \
        ${params.out_dir}/rsem_index/rsem_reference
    
    echo "RSEM index created"
    """
}