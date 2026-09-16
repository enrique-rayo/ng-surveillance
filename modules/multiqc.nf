process MULTIQC {
    tag "batch"
    container 'quay.io/biocontainers/multiqc:1.23--pyhdfd78af_0'
    publishDir 'results/multiqc', mode: 'copy'

    input:
    path fastqc_archives

    output:
    path "multiqc_report.html"
    path "multiqc_report_data"

    script:
    """
    multiqc --force --filename multiqc_report.html .
    """
}
