process SNIPPY {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/snippy:4.6.0--hdfd78af_2'
    publishDir 'results/snippy', mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    tuple path(reference), path(reference_indexes)

    output:
    tuple val(sample_id), path("${sample_id}_snippy"), emit: outdir
    tuple val(sample_id), path("${sample_id}_snippy/${sample_id}.raw.vcf.gz"), path("${sample_id}_snippy/${sample_id}.raw.vcf.gz.csi"), emit: raw_vcf

    script:
    def (forward, reverse) = reads
    """
    snippy \\
        --outdir ${sample_id}_snippy \\
        --ref ${reference} \\
        --R1 ${forward} \\
        --R2 ${reverse} \\
        --prefix ${sample_id} \\
        --mincov 5 \\
        --minfrac 0.75 \\
        --unmapped \\
        --cpus ${task.cpus}

    raw_vcf=${sample_id}_snippy/${sample_id}.raw.vcf
    compressed_vcf=\${raw_vcf}.gz

    if [[ -f "\${raw_vcf}" && ! -f "\${compressed_vcf}" ]]; then
        bgzip "\${raw_vcf}"
    fi

    if [[ ! -f "\${compressed_vcf}.csi" ]]; then
        bcftools index "\${compressed_vcf}"
    fi
    """
}

/*
 * Assembly-mapped branch (distinct from SNIPPY above, never merged with it):
 * each sample's reads are mapped to its own SPAdes assembly rather than a
 * shared external reference, with different, deliberately looser thresholds
 * (--mincov 10 --minfrac 0.05 vs. --mincov 5 --minfrac 0.75). This is a QC/
 * heterozygosity-detection branch, not a divergence-from-reference analysis:
 * expect very few filtered variants per sample (verified in practice: 6
 * filtered variants from 1279 raw candidate sites for one real sample),
 * with the raw candidate set feeding the het-site QC downstream.
 *
 * The reference here is per-sample (each sample's own clean contigs), unlike
 * SNIPPY's shared batch-wide reference + BWA index — hence a distinct process
 * rather than an aliased include, since the input shape itself differs, not
 * just the parameter values.
 */
process SNIPPY_ASSEMBLY {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/snippy:4.6.0--hdfd78af_2'
    publishDir 'results/snippy_assembly', mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(own_assembly)

    output:
    tuple val(sample_id), path("${sample_id}_snippy"), emit: outdir
    tuple val(sample_id), path("${sample_id}_snippy/${sample_id}.raw.vcf.gz"), path("${sample_id}_snippy/${sample_id}.raw.vcf.gz.csi"), emit: raw_vcf

    script:
    def (forward, reverse) = reads
    """
    snippy \\
        --outdir ${sample_id}_snippy \\
        --ref ${own_assembly} \\
        --R1 ${forward} \\
        --R2 ${reverse} \\
        --prefix ${sample_id} \\
        --mincov 10 \\
        --minfrac 0.05 \\
        --unmapped \\
        --cpus ${task.cpus}

    raw_vcf=${sample_id}_snippy/${sample_id}.raw.vcf
    compressed_vcf=\${raw_vcf}.gz

    if [[ -f "\${raw_vcf}" && ! -f "\${compressed_vcf}" ]]; then
        bgzip "\${raw_vcf}"
    fi

    if [[ ! -f "\${compressed_vcf}.csi" ]]; then
        bcftools index "\${compressed_vcf}"
    fi
    """
}

process SNIPPY_CORE {
    tag "$batch_name"
    // snippy-core builds one atomic alignment from every successful per-sample
    // directory. It has no option to skip an invalid input, so a failure must
    // terminate the run rather than publish a silently incomplete alignment.
    container 'quay.io/biocontainers/snippy:4.6.0--hdfd78af_2'
    publishDir { "results/snippy_core/${batch_name}" }, mode: 'copy'

    input:
    path snippy_dirs, arity: '1..*'
    path reference
    path mask
    val batch_name

    output:
    // snp-sites deliberately produces no .aln when the batch has zero SNPs.
    path "${batch_name}_snippy_masked.aln", emit: core_alignment, optional: true
    path "${batch_name}_snippy_masked.full.aln", emit: full_alignment
    path "${batch_name}_snippy_masked.vcf", emit: vcf
    path "${batch_name}_snippy_masked.tab", emit: tab
    path "${batch_name}_snippy_masked.txt", emit: stats
    path "${batch_name}_snippy_masked.ref.fa", emit: reference
    path "${batch_name}_snippy_masked.snippy-core.log", emit: log

    script:
    """
    set +e
    snippy-core \\
        --mask ${mask} \\
        --ref ${reference} \\
        --prefix ${batch_name}_snippy_masked \\
        ${snippy_dirs} \\
        2>&1 | tee ${batch_name}_snippy_masked.snippy-core.log
    snippy_core_status=\${PIPESTATUS[0]}
    set -e

    if [[ \${snippy_core_status} -ne 0 ]]; then
        if [[ \${snippy_core_status} -eq 2 ]] \\
            && grep -Fq 'Warning: No SNPs were detected so there is nothing to output.' ${batch_name}_snippy_masked.snippy-core.log \\
            && [[ -s ${batch_name}_snippy_masked.full.aln ]] \\
            && [[ -s ${batch_name}_snippy_masked.vcf ]] \\
            && ! grep -qv '^#' ${batch_name}_snippy_masked.vcf; then
            echo 'snippy-core produced no core-SNP alignment because this batch contains zero called variants.'
        else
            exit \${snippy_core_status}
        fi
    fi
    """
}
