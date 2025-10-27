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

process salmon_quant {
    tag "$SRR"
    
    publishDir "${params.outDir}/salmon_quant",
        mode: "copy",
        pattern: "*_quant.sf",
        enabled: params.pseudoAligner == "salmon"
    
    input:
    path salmon_index
    tuple val(SRR), path(fastq_files)

    output:
    path "${SRR}_quant.sf"
    tuple val(SRR), path("${SRR}_salmon_out/${SRR}_lib_format_counts.json"), emit: "lib_format_counts"
    tuple val(SRR), path("${SRR}_salmon_out/${SRR}_meta_info.json"), emit: "mapping_rate_info"

    module "Salmon"

    script:
    """
    salmon quant \\
        -i $salmon_index \\
        -l A \\
        -p $params.salmon_quant_threads \\
        -1 ${fastq_files[0]} \\
        -2 ${fastq_files[1]} \\
        -o ${SRR}_salmon_out
    mv ${SRR}_salmon_out/quant.sf ${SRR}_salmon_out/${SRR}_quant.sf
    mv ${SRR}_salmon_out/${SRR}_quant.sf ${SRR}_quant.sf
    mv ${SRR}_salmon_out/lib_format_counts.json ${SRR}_salmon_out/${SRR}_lib_format_counts.json
    mv ${SRR}_salmon_out/meta_info.json ${SRR}_salmon_out/${SRR}_meta_info.json
    """
}

process infer_strandness {
    tag "$SRR"

    input:
    tuple val(SRR), path(lib_format_counts)

    output:
    tuple val(SRR), path("${SRR}_strandness.txt")

    script:
    """
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
    " $lib_format_counts > ${SRR}_strandness.txt
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






process mapping_rate {
    tag "$SRR"

    input:
    tuple val(SRR), path(mapping_rate_info)

    output:
    tuple val(SRR), path("${SRR}_mapping_rate.txt")

    script:
    """
    jq -r "
        [
            .percent_mapped,
            .num_processed,
            .num_mapped,
            .num_decoy_fragments,
            .num_dovetail_fragments,
            .num_fragments_filtered_vm
        ] | @tsv
    " $mapping_rate_info > ${SRR}_mapping_rate.txt
    """
}

process merge_mapping_rates {
    tag "merge_mapping_rate_info"

    input:
    val srr_list
    path mapping_rate_files

    output:
    path "full_mapping_rate_info.tsv"

    script:
    """
    echo -e "Sample\\tmapping_rate\\tnum_processed_fragments\\tnum_mapped_fragments\\tnum_decoy_fragments\\tnum_dovetail_fragments\\tnum_fragments_filtered_vm" > full_mapping_rate_info.tsv

    mapping_files=(${mapping_rate_files.join(" ")})
    samples=(${srr_list.join(" ")})
    for idx in "\${!mapping_files[@]}"; do
        file="\${mapping_files[idx]}"
        sample="\${samples[idx]}"
        values=\$(cat "\$file")
        echo -e "\$sample\\t\$values" >> full_mapping_rate_info.tsv
    done
    """
}