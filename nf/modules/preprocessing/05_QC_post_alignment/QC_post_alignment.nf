nextflow.enable.dsl=2

process qc_post_alignment {
    // Define inputs:
    input:
    path star_sample_dir from file("${params.out_dir}/STAR_alignment/*/")
    path collapsed_gtf_file from params.rnaseQC.collapsed_gtf_file

    // Define output: rnaseQC quality control of the STAR alignment outputs
    output:
    path "${params.out_dir}/rnaseQC_metrics/${star_sample_dir.name}/*"

    // Define the script
    script:
    """
    # Create the log dir
    mkdir -p ${params.out_dir}/logs/qc_post_alignment
    
    # Load the module for RNAseQC
    module load RNA-SeQC

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample: ${star_sample_dir.name}"

    # Output directory
    out_dir="${params.out_dir}/rnaseQC_metrics/${star_sample_dir.name}"
    mkdir -p $out_dir

    # Star .bam file sorted by coordinates
    star_bam_file="${star_sample_dir}/${star_sample_dir.name}/${star_sample_dir.name}_Aligned.sortedByCoord.out.bam"

    rnaseqc ${collapsed_gtf_file} ${star_bam_file} --sample ${star_sample_dir.name} --stranded rf --verbose ${out_dir}/

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] RNAseQC of sample ${star_sample_dir.name} done"
    """
}