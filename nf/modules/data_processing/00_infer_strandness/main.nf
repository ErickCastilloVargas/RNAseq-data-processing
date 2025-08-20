#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process build_salmon_index {
    tag "salmon_index"

    input:
    path transcriptome

    output:
    path "salmon_index"

    module "Salmon"

    script:
    """
    salmon index -t $transcriptome -i salmon_index -k 31 --gencode
    """
}

process infer_strandness {
    tag "$SRR"

    input:
    path salmon_index
    tuple val(SRR), path(fastq_files)

    output:
    tuple val(SRR), path("${SRR}_strandness.txt")

    module "Salmon"

    script:
    """
    salmon quant \\
        -i salmon_index \\
        -l A \\
        -p $params.salmon_quant_threads \\
        -1 ${fastq_files[0]} \\
        -2 ${fastq_files[1]} \\
        -o ${SRR}_salmon_out

    jq -r "
        [
            .expected_format,
            .compatible_fragment_ratio,
            .num_compatible_fragments,
            .num_assigned_fragments,
            .num_frags_with_concordant_consistent_mappings,
            .num_frags_with_inconsistent_or_orphan_mappings,
            .strand_mapping_bias
        ] | @tsv
    " ${SRR}_salmon_out/lib_format_counts.json > ${SRR}_strandness.txt
    """
}

process merge_strandness_results {
    tag "merge_strandness_files"

    input:
    val srr_list
    path strandness_files

    output:
    path "full_strandness_inference.tsv"

    script:
    """
    echo -e "Sample\\texpected_format\\tcompatible_fragment_ratio\\tnum_compatible_fragments\\tnum_assigned_fragments\\tnum_consistent_mappings\\tnum_inconsistent_mappings\\tstrand_mapping_bias" > full_strandness_inference.tsv

    strand_files=(${strandness_files.join(" ")})
    samples=(${srr_list.join(" ")})
    for idx in "\${!strand_files[@]}"; do
        file="\${strand_files[idx]}"
        sample="\${samples[idx]}"
        values=\$(cat "\$file")
        echo -e "\$sample\\t\$values" >> full_strandness_inference.tsv
    done
    """
}