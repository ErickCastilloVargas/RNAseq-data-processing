nextflow.enable.dsl=2

process star_alignment {
    // Define input: List of directories (sample directories) within the trimmomatic output directory
    input:
    path sample_dir from file("${params.trimmomatic.out_dir}/*/")
    path star_index from params.index_building.out_dir_STAR_index
    val threads from params.STAR.threads

    // Define output: STAR alignment outputs
    output:
    path "${params.STAR.out_dir}/*"

    // Define the script
    script:
    """
    echo "Processing sample from directory: ${sample_dir}"

    # Output directory
    out_dir="${params.STAR.out_dir}/${sample_dir.name}"
    mkdir -p $out_dir

    # Input FASTQ files
    fastq1="${params.trimmomatic.out_dir}/${sample_dir.name}/${sample_dir.name}_1_paired.fastq.gz"
    fastq2="${params.trimmomatic.out_dir}/${sample_dir.name}/${sample_dir.name}_2_paired.fastq.gz"

    # Input genome index
    genome_index="$star_index

    # Do the STAR alignement 
    STAR --runMode alignReads \
        --genomeDir $genome_index \
        --runThreadN ${threads} \
        --twopassMode Basic \
        --genomeChrBinNbits 15 \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterType BySJout \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --limitSjdbInsertNsj 1200000 \
        --readFilesIn $fastq1 $fastq2 \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix ${out_dir}/${sample_dir.name}_ \
        --outSAMstrandField intronMotif \
        --alignSoftClipAtReferenceEnds Yes \
        --quantMode TranscriptomeSAM \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --genomeLoad NoSharedMemory \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType Junctions WithinBAM SoftClip \
        --chimMainSegmentMultNmax 1 \
        --outSAMattributes NH HI AS nM NM ch\
        --outSAMattrRGline ID:${sample_dir.name} SM:${sample_dir.name}

    echo "Sample ${sample_dir.name} aligned"
    """
}