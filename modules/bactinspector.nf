process BACTINSPECTOR {
    tag "$sample_id"
    errorStrategy 'ignore'
    container "${projectDir}/containers/bactinspector/bactinspector.sif"
    publishDir {
        task.process.contains('REASSEMBLY') \
            ? "results/bactinspector_reassembly/${batch_name}" \
            : "results/bactinspector/${batch_name}${task.process.contains('FILTERED') ? '/filtered' : ''}"
    }, mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)
    val batch_name

    output:
    tuple val(sample_id), path("${sample_id}_bactinspector.tsv"), emit: tsv

    script:
    """
    mkdir -p bactinspector_output
    ln -s ${contigs} ${sample_id}.fasta
    bactinspector check_species \
        -f ${sample_id}.fasta \
        -o bactinspector_output
    mv bactinspector_output/species_investigation_*.tsv \
        ${sample_id}_bactinspector.tsv
    """
}

process BACTINSPECTOR_BATCH_SUMMARY {
    tag "$batch_name"
    container "${projectDir}/containers/bactinspector/bactinspector.sif"
    publishDir {
        task.process.contains('REASSEMBLY') \
            ? "results/bactinspector_reassembly/${batch_name}" \
            : "results/bactinspector/${batch_name}${task.process.contains('FILTERED') ? '/filtered' : ''}"
    }, mode: 'copy'

    input:
    path sample_reports, arity: '0..*'
    val completed_sample_ids
    val expected_sample_ids
    val batch_name

    output:
    path "batch_summary.tsv", emit: tsv
    path "excluded_failed_samples.tsv", emit: excluded

    script:
    def shell_quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def expected = expected_sample_ids.collect(shell_quote).join(' ')
    def completed = completed_sample_ids.collect(shell_quote).join(' ')
    def concatenate_reports = sample_reports \
        ? "awk 'FNR == 1 && NR != 1 { next } { print }' ${sample_reports} > batch_summary.tsv" \
        : ': > batch_summary.tsv'
    """
    ${concatenate_reports}

    printf '%s\\n' ${expected} > expected_samples.txt
    printf '%s\\n' ${completed} > completed_samples.txt

    awk 'BEGIN {
             FS = OFS = "\\t"
             print "sample_id", "status", "detail"
         }
         NR == FNR { if (length(\$0)) completed[\$0] = 1; next }
         length(\$0) && !(\$0 in completed) {
             print \$0, "excluded_or_failed", \
                 "No BactInspector output was emitted; inspect the Nextflow run log and work directory for the upstream failure."
         }' completed_samples.txt expected_samples.txt \
         > excluded_failed_samples.tsv
    """
}
