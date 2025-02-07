nextflow.enable.dsl=2

process adapters_trimming {

    // Define input: List of directories (sample directories) within the main directory
    input:
    path sample_dir from file("${params.main_sample_dir}/*/")
    val threads from params.trimmomatic.threads
    path adapters_file from params.trimmomatic.adapters_file

    // Define output: Trimmomatic outputs
    output:
    path "${params.trimmomatic.out_dir}/${sample_dir.name}/*"

    // Define the script
    script:
    """
    # Load the module for Trimmomatic (this can be adjusted depending on your cluster setup)
    module load Trimmomatic

    echo "Processing sample from directory: ${sample_dir}"

    # Input FASTQ files (paired-end reads)
    fastq1="${sample_dir}/${sample_dir.name}_1.fastq.gz"
    fastq2="${sample_dir}/${sample_dir.name}_2.fastq.gz"
    
    # Output directory
    out_dir="${params.trimmomatic.out_dir}/${sample_dir.name}"
    mkdir -p $out_dir

    # Run Trimmomatic to remove adapters
    java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads ${threads} \
        $fastq1 \
        $fastq2 \
        $out_dir/${sample_dir.name}_1_paired.fastq.gz \
        $out_dir/${sample_dir.name}_1_unpaired.fastq.gz \
        $out_dir/${sample_dir.name}_2_paired.fastq.gz \
        $out_dir/${sample_dir.name}_2_unpaired.fastq.gz \
        ILLUMINACLIP:${adapters_file}:2:30:10 \
        HEADCROP:12 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    echo "All adapters are removed from ${sample_dir.name}!"
    """
}