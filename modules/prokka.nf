process PROKKA {
    tag "$sample_id"
    errorStrategy 'ignore'
    // Prokka probes the Java-based minced executable at startup. Concurrent
    // containers can contend for the same JVM perf-data file on this host.
    maxForks 1
    container 'quay.io/biocontainers/prokka:1.14.6--pl5321hdfd78af_5'
    publishDir 'results/prokka', mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)
    val genus

    output:
    tuple val(sample_id), path("${sample_id}_prokka"), emit: annotation
    tuple val(sample_id), path("${sample_id}_prokka/${sample_id}.gff"), emit: gff
    tuple val(sample_id), path("${sample_id}_prokka/${sample_id}.txt"), emit: summary

    script:
    def genus_arg = genus?.toString()?.trim() ? "--genus ${genus.toString().trim()}" : ''
    """
    prokka \
        --outdir ${sample_id}_prokka \
        --prefix ${sample_id} \
        --force \
        ${genus_arg} \
        ${contigs}
    """
}
