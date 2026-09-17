process ALIGNMENT_X_TO_N {
    tag "$batch_name"
    container 'quay.io/biocontainers/gubbins:3.3.5--py39pl5321he4a0461_0'
    publishDir { "results/phylogeny/${batch_name}/alignment_x_to_n" }, mode: 'copy'

    input:
    path full_alignment
    val batch_name

    output:
    path "${batch_name}_snippy_masked.full.x_to_n.aln", emit: alignment

    script:
    """
    tr 'X' 'N' < ${full_alignment} > ${batch_name}_snippy_masked.full.x_to_n.aln
    """
}

process GUBBINS {
    tag "$batch_name"
    // Gubbins produces one internally consistent, batch-wide recombination
    // analysis and tree. It cannot omit a failed isolate or publish a partial
    // batch result safely, so retain Nextflow's default terminate strategy.
    container 'quay.io/biocontainers/gubbins:3.3.5--py39pl5321he4a0461_0'
    publishDir { "results/phylogeny/${batch_name}/gubbins" }, mode: 'copy'

    input:
    path alignment
    val batch_name

    output:
    path "${batch_name}_snippy_masked.full.x_to_n.*", emit: results
    path "${batch_name}_snippy_masked.full.x_to_n.recombination_predictions.gff", emit: gff
    path "${batch_name}_snippy_masked.full.x_to_n.final_tree.tre", emit: tree

    script:
    """
    mkdir numba_cache_dir
    export NUMBA_CACHE_DIR="\$PWD/numba_cache_dir"

    run_gubbins.py \
        --threads ${task.cpus} \
        --first-tree-builder fasttree \
        --tree-builder iqtree \
        ${alignment}
    """
}

process MASK_GUBBINS_ALN {
    tag "$batch_name"
    // Like Gubbins and snippy-core, this is a single batch product. A missing
    // alignment/GFF or masking failure must terminate rather than look complete.
    container 'quay.io/biocontainers/gubbins:3.3.5--py39pl5321he4a0461_0'
    publishDir { "results/phylogeny/${batch_name}/masked_alignment" }, mode: 'copy'

    input:
    path core_alignment
    path recombination_gff
    val batch_name

    output:
    path "${batch_name}_snippy_masked.gubbins.aln", emit: alignment

    script:
    """
    mask_gubbins_aln.py \
        --aln ${core_alignment} \
        --gff ${recombination_gff} \
        --out ${batch_name}_snippy_masked.gubbins.aln
    """
}
