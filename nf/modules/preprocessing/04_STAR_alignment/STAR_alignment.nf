nextflow.enable.dsl=2

process star_alignment {
    // Define inputs:
    input:
    path sample_dir from file("${params.out_dir}/trimmomatic_trimming/*/")
    path star_index from ${params.out_dir}/star_index
    val threads from params.STAR.threads

    // Define output: STAR alignment outputs
    output:
    path "${params.out_dir}/STAR_alignment/${sample_dir.name}/*"

    // Define the script
    script:
    """
    # Create the log dir
    mkdir -p ${params.out_dir}/logs/star_alignment

    # Load the module for STAR
    module load STAR

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: ${sample_dir.name}"

    # Output directory
    out_dir="${params.out_dir}/STAR_alignment/${sample_dir.name}"
    mkdir -p $out_dir

    # Input FASTQ files
    fastq1="${sample_dir}/${sample_dir.name}/${sample_dir.name}_1_paired.fastq.gz"
    fastq2="${sample_dir}/${sample_dir.name}/${sample_dir.name}_2_paired.fastq.gz"

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

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Sample ${sample_dir.name} aligned"
    """
}