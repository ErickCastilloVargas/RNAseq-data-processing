nextflow.enable.dsl=2

process fastQC_post_trimming {
    // Define input: List of directories (sample directories) within the trimmomatic output directory
    input:
    path sample_dir from file("${params.trimmomatic.out_dir}/*/")
    val threads from params.fastQC.threads

    // Define output: FastQC outputs post trimming
    output:
    path "${params.fastQC.post_trimming_out_dir}/*"

    // Define the script
    script:
    """
    echo "Processing sample from directory: ${sample_dir}"

    # Output directory
    out_dir="${params.fastQC.post_trimming_out_dir}"
    mkdir -p $out_dir

    # Input FASTQ files
    fastq1="${params.trimmomatic.out_dir}/${sample_dir.name}/${sample_dir.name}_1_paired.fastq.gz"
    fastq2="${params.trimmomatic.out_dir}/${sample_dir.name}/${sample_dir.name}_2_paired.fastq.gz"

    fastqc -t ${threads} \
        ${fastq1} \
        ${fastq2} \
        -o ${out_dir}/

    echo "FastQC post trimming of sample ${sample_dir.name} done"
    """
}