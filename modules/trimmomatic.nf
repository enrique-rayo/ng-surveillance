process TRIMMOMATIC {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2'
    publishDir 'results/trimmomatic', mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_{forward,reverse}_paired.fq.gz"), emit: paired
    tuple val(sample_id), path("${sample_id}_{forward,reverse}_unpaired.fq.gz"), emit: unpaired

    script:
    """
    trimmomatic PE -threads ${task.cpus} \
        ${reads[0]} \
        ${reads[1]} \
        ${sample_id}_forward_paired.fq.gz \
        ${sample_id}_forward_unpaired.fq.gz \
        ${sample_id}_reverse_paired.fq.gz \
        ${sample_id}_reverse_unpaired.fq.gz \
        LEADING:25 \
        TRAILING:25 \
        MINLEN:36 \
        SLIDINGWINDOW:4:20
    """
}
