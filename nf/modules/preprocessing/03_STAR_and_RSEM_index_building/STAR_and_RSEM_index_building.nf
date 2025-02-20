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
    # Create the log dir 
    mkdir -p ${params.out_dir}/logs/star_and_RSEM_index_building

    # Load the module for STAR and RSEM
    module load STAR
    module load RSEM

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting to create STAR index"

    # STAR index generation
    STAR --runMode genomeGenerate \
        --genomeDir ${params.out_dir}/star_index \
        --genomeFastaFiles $reference_genome \
        --sjdbGTFfile $gtf_file \
        --sjdbOverhang 100 \
        --runThreadN $threads
     
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] STAR index created"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting to create RSEM index"

    # RSEM index generation
    rsem-prepare-reference --num-threads $threads \
        --gtf $gtf_file \
        $reference_genome \
        ${params.out_dir}/rsem_index/rsem_reference
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] RSEM index created"
    """
}