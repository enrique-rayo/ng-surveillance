process QUAST_BACTINSPECTOR_FLAGS {
    tag "$batch_name"
    container 'python:3.12.11-slim'
    publishDir { "results/quast_bactinspector_flags/${batch_name}" }, mode: 'copy'

    input:
    path quast_summary, name: 'quast_batch_summary.tsv'
    path quast_excluded, name: 'quast_excluded_failed_samples.tsv'
    path bactinspector_summary, name: 'bactinspector_batch_summary.tsv'
    path bactinspector_excluded, name: 'bactinspector_excluded_failed_samples.tsv'
    val expected_sample_ids
    val batch_name
    path flags_script
    path qc_flagging_script

    output:
    path "flagged_sample_ids.txt", emit: flagged_ids

    script:
    def expected = expected_sample_ids.join(',')
    """
    python ${flags_script} \
        --quast-summary quast_batch_summary.tsv \
        --quast-excluded quast_excluded_failed_samples.tsv \
        --bactinspector-summary bactinspector_batch_summary.tsv \
        --bactinspector-excluded bactinspector_excluded_failed_samples.tsv \
        --expected-samples '${expected}' \
        > flagged_sample_ids.txt
    """
}
