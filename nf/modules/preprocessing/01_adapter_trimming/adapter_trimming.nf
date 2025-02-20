nextflow.enable.dsl=2

process adapters_trimming {
    // Define inputs:
    input:
    path sample_dir from file("${params.main_sample_dir}/*/")
    val threads from params.trimmomatic.threads
    path adapters_file from params.trimmomatic.adapters_file

    // Define output: adapters and low quality read fragments removal from samples 
    output:
    path "${params.out_dir}/trimmomatic_trimming/${sample_dir.name}/*"

    // Define the script
    script:
    """
    # Create the log dir
    mkdir -p ${params.out_dir}/logs/adapters_trimming

    # Load the module for Trimmomatic
    module load Trimmomatic

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: ${sample_dir.name}"

    # Input FASTQ files (paired-end reads)
    fastq1="${sample_dir}/${sample_dir.name}_1.fastq.gz"
    fastq2="${sample_dir}/${sample_dir.name}_2.fastq.gz"
    
    # Output directory
    out_dir="${params.out_dir}/trimmomatic_trimming/${sample_dir.name}"
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
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] All adapters are removed from ${sample_dir.name}!"
    """
}