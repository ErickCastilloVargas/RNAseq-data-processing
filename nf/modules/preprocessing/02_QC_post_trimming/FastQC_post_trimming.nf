nextflow.enable.dsl=2

process fastQC_post_trimming {
    // Define inputs:
    input:
    path sample_dir from file("${params.out_dir}/trimmomatic_trimming/*/")
    val threads from params.fastQC.threads

    // Define output: FastQC outputs post trimming
    output:
    path "${params.out_dir}/fastQC_metrics/${sample_dir.name}/*"

    // Define the script
    script:
    """
    # Create the log dir 
    mkdir -p ${params.out_dir}/logs/fastQC_post_trimming
    
    # Load the module for fastQC
    module load FastQC

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: ${sample_dir.name}"

    # Output directory
    out_dir="${params.out_dir}/fastQC_metrics/${sample_dir.name}"
    mkdir -p $out_dir

    # Input FASTQ files
    fastq1="${params.trimmomatic.out_dir}/${sample_dir.name}/${sample_dir.name}_1_paired.fastq.gz"
    fastq2="${params.trimmomatic.out_dir}/${sample_dir.name}/${sample_dir.name}_2_paired.fastq.gz"

    fastqc -t ${threads} \
        ${fastq1} \
        ${fastq2} \
        -o ${out_dir}/

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] FastQC post trimming of sample ${sample_dir.name} done"
    """
}