process QUAST {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/quast:5.3.0--py311pl5321hc84137b_1'
    publishDir {
        task.process.contains('REASSEMBLY') \
            ? "results/quast_reassembly/${batch_name}" \
            : "results/quast/${batch_name}${task.process.contains('FILTERED') ? '/filtered' : ''}"
    }, mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)
    val batch_name

    output:
    tuple val(sample_id), path("${sample_id}"), emit: report
    tuple val(sample_id), path("${sample_id}_transposed_report.tsv"), emit: tsv

    script:
    """
    ln -s ${contigs} ${sample_id}.fasta
    quast.py ${sample_id}.fasta \
        -o ${sample_id} \
        --threads ${task.cpus}
    ln -s ${sample_id}/transposed_report.tsv \
        ${sample_id}_transposed_report.tsv
    """
}

process QUAST_BATCH_SUMMARY {
    tag "$batch_name"
    container 'quay.io/biocontainers/quast:5.3.0--py311pl5321hc84137b_1'
    publishDir {
        task.process.contains('REASSEMBLY') \
            ? "results/quast_reassembly/${batch_name}" \
            : "results/quast/${batch_name}${task.process.contains('FILTERED') ? '/filtered' : ''}"
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
                 "No QUAST output was emitted; inspect the Nextflow run log and work directory for the upstream failure."
         }' completed_samples.txt expected_samples.txt \
         > excluded_failed_samples.tsv
    """
}
