process ASSEMBLY_SELECTION {
    tag "$batch_name"
    container 'debian:12.11-slim'
    publishDir { "results/assembly_selection/${batch_name}" }, mode: 'copy'

    input:
    val selection_rows
    val batch_name

    output:
    path "assembly_selection.tsv", emit: manifest

    script:
    def shell_quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def rows = selection_rows
        .sort { left, right -> left[0].toString() <=> right[0].toString() }
        .collect { row -> row.collect { it.toString() }.join('\t') }
        .collect(shell_quote)
        .join(' ')
    def append_rows = rows \
        ? "printf '%s\\n' ${rows} >> assembly_selection.tsv" \
        : ':'
    """
    printf 'sample_id\tselected_source\tpromoted_for_prokka_pyngost\n' \
        > assembly_selection.tsv
    ${append_rows}
    """
}
