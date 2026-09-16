/*
 * Filter both paired reads with the target-taxon read-ID list using
 * seqtk subseq. Consumes Trimmomatic's paired outputs and gzips both
 * results to match the pipeline's FASTQ convention.
 */
process SEQTK_NG_FILTER {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/seqtk:1.5--h577a1d6_1'
    publishDir 'results/seqtk_ng_filter', mode: 'copy'

    input:
    tuple val(sample_id), path(forward_paired), path(reverse_paired), path(ng_list)

    output:
    tuple val(sample_id), path("${sample_id}_NG_gen_R1.fq.gz"), path("${sample_id}_NG_gen_R2.fq.gz"), emit: ng_filtered_reads

    script:
    """
    awk '{ print; print \$0 "/1" }' ${ng_list} > forward_ng_ids.lst
    awk '{ print; print \$0 "/2" }' ${ng_list} > reverse_ng_ids.lst
    seqtk subseq ${forward_paired} forward_ng_ids.lst > ${sample_id}_NG_gen_R1.fq
    seqtk subseq ${reverse_paired} reverse_ng_ids.lst > ${sample_id}_NG_gen_R2.fq
    gzip ${sample_id}_NG_gen_R1.fq ${sample_id}_NG_gen_R2.fq
    """
}
