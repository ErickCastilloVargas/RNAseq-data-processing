#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process star_alignment {
    clusterOptions = { 
        "--output=STAR_${SRR}.out --error=STAR_${SRR}.err" 
    }
    publishDir "results/STAR_alignment", pattern: "*.out.bam"
    publishDir "results/logs/STAR", pattern: "*.{out,err}"

    input:
    tuple val(SRR), path(fastq_files)
    path star_index

    output:
    tuple val(SRR), path("*.out.bam"), emit: "bam_files"
    path "*.{out,err}"

    script:
    """
    module load STAR

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: ${SRR}"

    STAR --runMode alignReads \\
        --genomeDir $star_index \\
        --runThreadN ${params.STAR.threads} \\
        --twopassMode Basic \\
        --outFilterMultimapNmax 20 \\
        --alignSJoverhangMin 8 \\
        --alignSJDBoverhangMin 1 \\
        --outFilterMismatchNmax 999 \\
        --outFilterMismatchNoverLmax 0.1 \\
        --alignIntronMin 20 \\
        --alignIntronMax 1000000 \\
        --alignMatesGapMax 1000000 \\
        --outFilterType BySJout \\
        --outFilterScoreMinOverLread 0.33 \\
        --outFilterMatchNminOverLread 0.33 \\
        --limitSjdbInsertNsj 1200000 \\
        --readFilesIn $fastq_files \\
        --readFilesCommand gunzip -c \\
        --outFileNamePrefix ${SRR}_ \\
        --outSAMstrandField intronMotif \\
        --alignSoftClipAtReferenceEnds Yes \\
        --quantMode TranscriptomeSAM \\
        --outSAMtype BAM Unsorted \\
        --outSAMunmapped Within \\
        --genomeLoad NoSharedMemory \\
        --chimSegmentMin 15 \\
        --chimJunctionOverhangMin 15 \\
        --chimOutType Junctions WithinBAM SoftClip \\
        --chimMainSegmentMultNmax 1 \\
        --outSAMattributes NH HI AS nM NM ch \\
        --outSAMattrRGline ID:${SRR} SM:${SRR}

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Sample ${SRR} aligned"
    """
}

process star_alignment_no_adapter_trimming {
    clusterOptions = { 
        "--output=STAR_${SRR}.out --error=STAR_${SRR}.err" 
    }
    publishDir "results/STAR_alignment", pattern: "*.out.bam"
    publishDir "results/logs/STAR", pattern: "*.{out,err}"

    // Define inputs:
    input:
    tuple val(SRR), path(fastq_files)
    path star_index

    output:
    tuple val(SRR), path("*.out.bam"), emit: "bam_files"
    path "*.{out,err}"

    // Define the script
    script:
    """
    module load STAR

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: ${SRR}"

    STAR --runMode alignReads \\
        --genomeDir $star_index \\
        --runThreadN ${params.STAR.threads} \\
        --twopassMode Basic \\
        --outFilterMultimapNmax 20 \\
        --alignSJoverhangMin 8 \\
        --alignSJDBoverhangMin 1 \\
        --outFilterMismatchNmax 999 \\
        --outFilterMismatchNoverLmax 0.1 \\
        --alignIntronMin 20 \\
        --alignIntronMax 1000000 \\
        --alignMatesGapMax 1000000 \\
        --outFilterType BySJout \\
        --outFilterScoreMinOverLread 0.33 \\
        --outFilterMatchNminOverLread 0.33 \\
        --limitSjdbInsertNsj 1200000 \\
        --readFilesIn $fastq_files \\
        --readFilesCommand gunzip -c \\
        --outFileNamePrefix ${SRR}_ \\
        --outSAMstrandField intronMotif \\
        --alignSoftClipAtReferenceEnds Yes \\
        --quantMode TranscriptomeSAM \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --genomeLoad NoSharedMemory \\
        --chimSegmentMin 15 \\
        --chimJunctionOverhangMin 15 \\
        --chimOutType Junctions WithinBAM SoftClip \\
        --chimMainSegmentMultNmax 1 \\
        --outSAMattributes NH HI AS nM NM ch \\
        --outSAMattrRGline ID:${SRR} SM:${SRR}

    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Sample ${SRR} aligned"
    """
}