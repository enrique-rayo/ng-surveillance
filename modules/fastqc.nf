process FASTQC {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'biocontainers/fastqc:v0.11.9_cv8'
    publishDir 'results/fastqc', mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), emit: zip
    tuple val(sample_id), path("*_fastqc.html"), emit: html

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}
