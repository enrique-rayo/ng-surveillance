process KRAKEN2 {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_1'
    publishDir 'results/kraken2', mode: 'copy'
    maxForks 1
    memory '90 GB'

    input:
    tuple val(sample_id), path(reads)
    path kraken2_db

    output:
    tuple val(sample_id), path("${sample_id}_kraken_report.txt"), emit: report
    tuple val(sample_id), path("${sample_id}_kraken_results.txt"), emit: results
    tuple val(sample_id), path("${sample_id}_classified{_1,_2}.fq"), emit: classified
    tuple val(sample_id), path("${sample_id}_unclassified{_1,_2}.fq"), emit: unclassified

    script:
    """
    kraken2 --db ${kraken2_db} \
        --threads ${task.cpus} \
        --paired \
        --gzip-compressed \
        --report-zero-counts \
        --output ${sample_id}_kraken_results.txt \
        --report ${sample_id}_kraken_report.txt \
        --classified-out '${sample_id}_classified#.fq' \
        --unclassified-out '${sample_id}_unclassified#.fq' \
        ${reads[0]} \
        ${reads[1]}
    """
}

process KRAKEN2_TAXID_SUMMARY {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'python:3.12.11-slim'
    publishDir 'results/kraken2', mode: 'copy'

    input:
    tuple val(sample_id), path(report)
    path summary_script
    val target_genus_taxid
    val target_species_taxid

    output:
    tuple val(sample_id), path("${sample_id}_kraken_summary.csv"), emit: csv

    script:
    """
    python ${summary_script} \
        --report ${report} \
        --sample-id ${sample_id} \
        --output ${sample_id}_kraken_summary.csv \
        --target-genus-taxid ${target_genus_taxid} \
        --target-species-taxid ${target_species_taxid}
    """
}

process KRAKEN2_BATCH_SUMMARY {
    tag "$batch_name"
    container 'quay.io/biocontainers/matplotlib:3.1.2--2'
    publishDir 'results/kraken2/kraken_summary', mode: 'copy'

    input:
    path reports
    path taxid_summaries
    val batch_name
    path batch_summary_script
    path qc_flagging_script

    output:
    path "kraken2_taxa_summary.csv", emit: csv
    path "kraken2_batch_summary.html", emit: html

    script:
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    python ${batch_summary_script} \
        --results-dir . \
        --batch-name '${batch_name}' \
        --csv-output kraken2_taxa_summary.csv \
        --html-output kraken2_batch_summary.html
    """
}
