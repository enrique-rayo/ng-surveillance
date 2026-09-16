/*
 * Adapted from nf-core/modules ariba/run.
 *
 * The upstream module expects a prepared ARIBA database as a .tar.gz archive.
 * Extract into an explicit directory so both flat archives and archives with a
 * top-level directory have a stable path at runtime.
 */
process ARIBA_RUN {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/ariba:2.14.6--py39h67e14b5_3'
    publishDir 'results/ariba', mode: 'copy'
    // ARIBA 2.14.6 uses a multiprocessing manager socket even with --threads 1;
    // concurrent samples can collide on its listening address.
    maxForks 1

    input:
    tuple val(sample_id), path(reads)
    path db

    output:
    tuple val(sample_id), path("${sample_id}"), emit: results

    script:
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    mkdir ariba_db
    tar -xzf ${db} -C ariba_db
    ariba run \
        --verbose \
        --threads ${task.cpus} \
        ariba_db \
        ${reads[0]} \
        ${reads[1]} \
        ${sample_id}
    """
}

/*
 * Aggregate every successful per-sample ARIBA run with the exact command used
 * by SensiTyper. Keep a self-contained reports/ tree beside the summary so
 * the CSV's first-column paths remain valid when both outputs are staged into
 * a downstream Nextflow task or copied to results/ariba_summary.
 *
 * NG-specific by design: the two sed rewrites below correct known ARIBA
 * report artefacts for two specific penA mosaic-allele insertion variants
 * (D147_T148insT, R146_D147insR). These are validated scientific corrections
 * for N. gonorrhoeae penA typing, not generic tunables — they will not apply,
 * and should not be adapted, for another organism's ARIBA output.
 */
process ARIBA_SUMMARY {
    tag "batch"
    container 'quay.io/biocontainers/ariba:2.14.6--py39h67e14b5_3'
    publishDir 'results/ariba_summary', mode: 'copy'

    input:
    path ariba_results

    output:
    path "ariba_summary.csv", emit: csv
    path "reports", emit: reports

    script:
    """
    mkdir -p reports
    : > filenames.txt

    for result_dir in ${ariba_results}; do
        sample_id=\$(basename "\$result_dir")
        test -r "\$result_dir/report.tsv"
        report_dir="reports/\${sample_id}_ARIBA"
        mkdir -p "\$report_dir"

        # SensiTyper's ariba_batch_v0.2.py creates report_complete.tsv with
        # these two exact penA insertion rewrites before running the summary.
        sed \
            -e '/D147_T148insT/ { s/0\\t\\.\\tp\\t\\.\\t0\\tD147_T148insT/1\\tSNP\\tp\\tD147_T148insT\\t1\\tD147_T148insT/; b; }' \
            -e '/R146_D147insR/ s/0\\t\\.\\tp\\t\\.\\t0\\tR146_D147insR/1\\tSNP\\tp\\tR146_D147insR\\t1\\tR146_D147insR/' \
            "\$result_dir/report.tsv" > "\$report_dir/report_complete.tsv"

        printf 'reports/%s_ARIBA/report_complete.tsv\\n' "\$sample_id" >> filenames.txt
    done

    test -s filenames.txt
    sort -o filenames.txt filenames.txt

    export MPLCONFIGDIR="\$PWD/.matplotlib"
    ariba summary ariba_summary \
        -f filenames.txt \
        --cluster_cols assembled,ref_seq,pct_id \
        --col_filter n \
        --row_filter n \
        --no_tree \
        --v_groups \
        --known_variants
    """
}
