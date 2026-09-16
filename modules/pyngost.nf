/*
 * NG-specific by design: pyngoST's -s NG-STAR,MLST flag requests typing
 * schemes defined specifically for N. gonorrhoeae (NG-STAR = N. gonorrhoeae
 * Sequence Typing for Antimicrobial Resistance). Adapting this pipeline for
 * another organism means replacing this module with that organism's own
 * typing tool, not adjusting the scheme flag.
 */
process PYNGOST {
    tag "$batch_name"
    // This is one atomic batch report: a process failure must fail the run rather
    // than silently omit typing for every sample, so retain Nextflow's default
    // terminate error strategy.
    container 'quay.io/biocontainers/pyngost:1.1.3--pyh7e72e81_0'
    publishDir { "results/pyngost/${batch_name}" }, mode: 'copy'

    input:
    path assemblies, arity: '1..*'
    path db
    val batch_name
    val excluded_sample_ids

    output:
    path "pyngost_results", emit: results

    script:
    def exclusion_manifest = excluded_sample_ids
        ? """cat > pyngost_results/excluded_from_pyngost.tsv <<'EOF'
sample_id\treason
${excluded_sample_ids.sort().collect { "${it}\tmanually excluded via params.pyngost_exclude_samples — see docs/OPEN_PARAMETERS.md" }.join('\n')}
EOF
"""
        : ''

    """
    mkdir pyngost_results
    ${exclusion_manifest}
    pyngoST.py \
        -i ${assemblies} \
        -p ${db} \
        -s NG-STAR,MLST \
        -c \
        -m \
        -b \
        -a \
        -q pyngost_results \
        -o ${batch_name}_pyngoST_out.txt
    """
}
