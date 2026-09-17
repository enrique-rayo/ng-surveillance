include { SNIPPY; SNIPPY_ASSEMBLY; SNIPPY_CORE } from '../modules/snippy'
include { HET_SITES_FILTER; HET_SITES_PANEL; HET_SITES_QUERY; HET_SITES_BIN; HET_SITES_BATCH_SUMMARY } from '../modules/het_sites'
include { ALIGNMENT_X_TO_N; GUBBINS; MASK_GUBBINS_ALN } from '../modules/gubbins'

workflow PHYLOGENY {
    take:
    trimmed_paired_ch            // tuple(sample_id, [forward_paired, reverse_paired])
    assembly_for_downstream_ch   // tuple(sample_id, clean_contigs)
    expected_sample_ids_ch       // list of all sample IDs in the batch
    batch_name                   // val

    main:
    qc_control_pattern = java.util.regex.Pattern.compile(
        params.qc_control_sample_pattern ?: '(?i)^NGIVPZEC\\d+(?:_S\\d+)?$'
    )
    run_phylogeny = (params.run_phylogeny ?: false).toString().toBoolean()
    run_assembly_mapped_snippy = (params.run_assembly_mapped_snippy ?: false).toString().toBoolean()

    // --- Reference-mapped branch (FA1090) ---
    if (params.snippy_reference) {
        snippy_reference = file(params.snippy_reference, checkIfExists: true)
        snippy_reference_indexes = ['.amb', '.ann', '.bwt', '.pac', '.sa']
            .collect { suffix -> file("${params.snippy_reference}${suffix}", checkIfExists: true) }
        snippy_reference_ch = Channel.value([snippy_reference, snippy_reference_indexes])
        SNIPPY(trimmed_paired_ch, snippy_reference_ch)

        HET_SITES_FILTER(SNIPPY.out.raw_vcf, batch_name)
        HET_SITES_PANEL(
            HET_SITES_FILTER.out.sites.map { sample_id, sites -> sites }.collect(),
            batch_name
        )

        het_successful_raw_vcf_ch = SNIPPY.out.raw_vcf
            .join(HET_SITES_FILTER.out.sites)
            .map { sample_id, raw_vcf, raw_vcf_index, sites ->
                tuple(sample_id, raw_vcf, raw_vcf_index)
            }

        HET_SITES_QUERY(
            het_successful_raw_vcf_ch,
            HET_SITES_PANEL.out.panel,
            batch_name
        )
        HET_SITES_BIN(
            HET_SITES_QUERY.out.queried,
            Channel.value(file("${projectDir}/bin/het_per_100-mer.py", checkIfExists: true)),
            batch_name
        )
        HET_SITES_BATCH_SUMMARY(
            HET_SITES_BIN.out.summary.map { sample_id, summary -> summary }.collect().ifEmpty([]),
            HET_SITES_BIN.out.summary.map { sample_id, summary -> sample_id }.collect().ifEmpty([]),
            expected_sample_ids_ch,
            batch_name
        )

        if (params.snippy_mask) {
            snippy_core_exclude_sample_ids = (params.snippy_core_exclude_samples ?: '')
                .toString()
                .split(',')
                .collect { it.trim() }
                .findAll { it }
                .toSet()

            snippy_core_outdirs_ch = SNIPPY.out.outdir
                .filter { sample_id, snippy_dir ->
                    !qc_control_pattern.matcher(sample_id).matches() &&
                        !snippy_core_exclude_sample_ids.contains(sample_id)
                }

            SNIPPY_CORE(
                snippy_core_outdirs_ch.map { sample_id, snippy_dir -> snippy_dir }.collect(),
                Channel.value(snippy_reference),
                Channel.value(file(params.snippy_mask, checkIfExists: true)),
                batch_name
            )

            if (run_phylogeny) {
                ALIGNMENT_X_TO_N(SNIPPY_CORE.out.full_alignment, batch_name)
                GUBBINS(ALIGNMENT_X_TO_N.out.alignment, batch_name)
                MASK_GUBBINS_ALN(
                    SNIPPY_CORE.out.core_alignment,
                    GUBBINS.out.gff,
                    batch_name
                )
            }
        } else {
            log.warn 'Snippy-core skipped: provide the required --snippy_mask parameter (BED file) to build the batch alignment.'
            if (run_phylogeny) {
                log.warn 'Phylogeny skipped: --run_phylogeny requires --snippy_mask so snippy-core can build its input alignments.'
            }
        }
    } else {
        log.warn 'Reference-mapped Snippy skipped: provide the required --snippy_reference parameter (FASTA with adjacent BWA index files).'
        if (run_phylogeny) {
            log.warn 'Phylogeny skipped: --run_phylogeny requires --snippy_reference and --snippy_mask.'
        }
    }

    // --- Assembly-mapped branch (each sample vs. its own assembly) ---
    // Opt-in: unlike the reference-mapped branch above, this needs no external
    // database, so there is nothing to naturally gate it on. Its purpose is
    // heterozygous/mixed-site QC, not divergence analysis — expect very few
    // filtered variants per sample (validated on one real sample: 6 filtered
    // from 1279 raw candidates). Het-site QC for this branch (filter/panel/
    // query/bin) is NOT wired here — deliberately deferred, since
    // het_per_100-mer.py's fixed-bin approach is specific to the single,
    // fixed-length FA1090 reference and does not apply to each sample's own,
    // differently-sized, multi-contig assembly. See modules/het_sites.nf for
    // detail. A per-sample-assembly-appropriate binning/QC method needs to be
    // designed before this branch's het-site QC can be added.
    if (run_assembly_mapped_snippy) {
        snippy_assembly_input_ch = trimmed_paired_ch.join(assembly_for_downstream_ch)
        SNIPPY_ASSEMBLY(snippy_assembly_input_ch)
    } else {
        log.warn 'Assembly-mapped Snippy skipped: set --run_assembly_mapped_snippy to run it (het-site QC for this branch is not yet implemented).'
    }
}
