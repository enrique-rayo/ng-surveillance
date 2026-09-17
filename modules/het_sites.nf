/*
 * Reference-mapped branch only. HET_SITES_BIN below relies on
 * het_per_100-mer.py's fixed 100-bin structure sized specifically to the
 * FA1090 reference genome's length, and discards VCF CHROM entirely
 * (assumes one single-chromosome reference). Both assumptions hold for the
 * reference-mapped branch (every sample shares the same FA1090 reference)
 * but break for the assembly-mapped branch, where each sample has its own
 * assembly of a different length, spread across multiple contigs.
 *
 * HET_SITES_FILTER/PANEL/QUERY below are otherwise generic (pure VCF/
 * bcftools operations, no reference-specific assumptions) and could in
 * principle run against assembly-mapped VCFs too — but the assembly-mapped
 * branch's het-site QC (a binning/summary methodology appropriate to
 * per-sample, multi-contig assemblies) has not yet been designed and is
 * deliberately deferred rather than invented here. Extending this to the
 * assembly-mapped branch needs a project-lead decision on binning approach
 * before implementation, not a silent reuse of this fixed-bin script.
 */
process HET_SITES_FILTER {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/snippy:4.6.0--hdfd78af_2'
    publishDir { "results/het_sites/${batch_name}/per_sample" }, mode: 'copy'

    input:
    tuple val(sample_id), path(raw_vcf), path(raw_vcf_index)
    val batch_name

    output:
    tuple val(sample_id), path("${sample_id}_hets.bcf"), emit: bcf
    tuple val(sample_id), path("${sample_id}_het_sites.txt"), emit: sites

    script:
    """
    bcftools annotate \\
        -x '^INFO/DP' \\
        -i 'INFO/AF>0.1 && INFO/AF<0.9 && INFO/DP>10' \\
        ${raw_vcf} \\
        -o ${sample_id}_hets.bcf \\
        -O b

    bcftools query \\
        -f '%CHROM\t%POS\n' \\
        ${sample_id}_hets.bcf \\
        > ${sample_id}_het_sites.txt
    """
}

process HET_SITES_PANEL {
    tag "$batch_name"
    container 'quay.io/biocontainers/snippy:4.6.0--hdfd78af_2'
    publishDir { "results/het_sites/${batch_name}" }, mode: 'copy'

    input:
    path sample_sites, arity: '0..*'
    val batch_name

    output:
    path "${batch_name}_unique_het_sites.txt", emit: panel

    script:
    def build_panel = sample_sites \
        ? "cat ${sample_sites} | sort -V | uniq > ${batch_name}_unique_het_sites.txt" \
        : ": > ${batch_name}_unique_het_sites.txt"
    """
    ${build_panel}
    """
}

process HET_SITES_QUERY {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/snippy:4.6.0--hdfd78af_2'
    publishDir { "results/het_sites/${batch_name}/per_sample" }, mode: 'copy'

    input:
    tuple val(sample_id), path(raw_vcf), path(raw_vcf_index)
    path panel
    val batch_name

    output:
    tuple val(sample_id), path("${sample_id}_all_het_sites.txt"), emit: queried

    script:
    """
    if [[ -s ${panel} ]]; then
        bcftools query \\
            -R ${panel} \\
            -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n' \\
            ${raw_vcf} \\
            -o ${sample_id}_all_het_sites.txt
    else
        : > ${sample_id}_all_het_sites.txt
    fi
    """
}

process HET_SITES_BIN {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'python:3.12.11-slim'
    publishDir { "results/het_sites/${batch_name}/per_sample" }, mode: 'copy'

    input:
    tuple val(sample_id), path(queried_sites)
    path bin_script
    val batch_name

    output:
    tuple val(sample_id), path("${sample_id}_het_bin_summary.tsv"), emit: summary

    script:
    """
    python ${bin_script} \\
        ${queried_sites} \\
        --sample-id ${sample_id} \\
        --output ${sample_id}_het_bin_summary.tsv
    """
}

process HET_SITES_BATCH_SUMMARY {
    tag "$batch_name"
    container 'python:3.12.11-slim'
    publishDir { "results/het_sites/${batch_name}" }, mode: 'copy'

    input:
    path sample_summaries, arity: '0..*'
    val completed_sample_ids
    val expected_sample_ids
    val batch_name

    output:
    path "${batch_name}_het_sites_summary.tsv", emit: summary
    path "excluded_failed_samples.tsv", emit: excluded

    script:
    def shell_quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def expected = expected_sample_ids.collect(shell_quote).join(' ')
    def completed = completed_sample_ids.collect(shell_quote).join(' ')
    def append_summaries = sample_summaries \
        ? "cat ${sample_summaries} | sort -k1,1 -k2,2n >> ${batch_name}_het_sites_summary.tsv" \
        : ':'
    """
    printf 'Sample\tNum. het. sites\tNum. bins\n' > ${batch_name}_het_sites_summary.tsv
    ${append_summaries}

    printf '%s\n' ${expected} > expected_samples.txt
    printf '%s\n' ${completed} > completed_samples.txt

    awk 'BEGIN {
             FS = OFS = "\t"
             print "sample_id", "status", "detail"
         }
         NR == FNR { if (length(\$0)) completed[\$0] = 1; next }
         length(\$0) && !(\$0 in completed) {
             print \$0, "excluded_or_failed", \
                 "No final het-site bin summary was emitted; inspect the Nextflow log and upstream HET_SITES work directories."
         }' completed_samples.txt expected_samples.txt \
         > excluded_failed_samples.tsv
    """
}
