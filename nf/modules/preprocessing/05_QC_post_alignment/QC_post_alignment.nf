nextflow.enable.dsl=2

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
    # Load the module for RNAseQC
    module load RNA-SeQC

    echo "Processing sample from directory: ${star_sample_dir}"

    # Output directory
    out_dir="${params.out_dir}/rnaseQC_metrics/${star_sample_dir.name}"
    mkdir -p $out_dir

    # Star .bam file sorted by coordinates
    star_bam_file="${star_sample_dir}/${star_sample_dir.name}/${star_sample_dir.name}_Aligned.sortedByCoord.out.bam"

    rnaseqc ${collapsed_gtf_file} ${star_bam_file} --sample ${star_sample_dir.name} --stranded rf --verbose ${out_dir}/

    echo "RNAseQC of sample ${star_sample_dir.name} done"
    """