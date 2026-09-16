process REDUCE_CONTIGS {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/seqtk:1.5--h577a1d6_1'
    publishDir {
        task.process.contains('REASSEMBLY') \
            ? 'results/spades_reassembly/clean_contigs' \
            : 'results/spades/clean_contigs'
    }, mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)
    val min_contig_length

    output:
    tuple val(sample_id), path("${sample_id}_clean.fasta"), emit: clean_contigs

    script:
    """
    seqtk seq -L ${min_contig_length} ${contigs} > ${sample_id}_clean.fasta
    test -s ${sample_id}_clean.fasta
    """
}
