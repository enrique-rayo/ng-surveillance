include { FASTQC } from '../modules/fastqc'
include { MULTIQC } from '../modules/multiqc'
include { TRIMMOMATIC } from '../modules/trimmomatic'
include { KRAKEN2; KRAKEN2_TAXID_SUMMARY; KRAKEN2_BATCH_SUMMARY } from '../modules/kraken2'
include { SPADES } from '../modules/spades'
include { SPADES as SPADES_REASSEMBLY } from '../modules/spades'
include { REDUCE_CONTIGS } from '../modules/reduce_contigs'
include { REDUCE_CONTIGS as REDUCE_CONTIGS_REASSEMBLY } from '../modules/reduce_contigs'
include { BACTINSPECTOR; BACTINSPECTOR_BATCH_SUMMARY } from '../modules/bactinspector'
include { BACTINSPECTOR as BACTINSPECTOR_FILTERED; BACTINSPECTOR_BATCH_SUMMARY as BACTINSPECTOR_FILTERED_BATCH_SUMMARY } from '../modules/bactinspector'
include { BACTINSPECTOR as BACTINSPECTOR_REASSEMBLY } from '../modules/bactinspector'
include { QUAST; QUAST_BATCH_SUMMARY } from '../modules/quast'
include { QUAST as QUAST_FILTERED; QUAST_BATCH_SUMMARY as QUAST_FILTERED_BATCH_SUMMARY } from '../modules/quast'
include { QUAST as QUAST_REASSEMBLY } from '../modules/quast'
include { QUAST_BACTINSPECTOR_FLAGS } from '../modules/quast_bactinspector_flags'
include { NG_READ_LIST } from '../modules/ng_read_list'
include { SEQTK_NG_FILTER } from '../modules/seqtk_ng_filter'
include { ASSEMBLY_SELECTION } from '../modules/assembly_selection'

workflow QC_ASSEMBLY {
    take:
    read_pairs_ch   // tuple(sample_id, [R1, R2]) from Channel.fromFilePairs
    batch_name      // val

    main:
    expected_sample_ids_ch = read_pairs_ch
        .map { sample_id, reads -> sample_id }
        .collect()

    FASTQC(read_pairs_ch)
    MULTIQC(FASTQC.out.zip.map { sample_id, archive -> archive }.collect())
    TRIMMOMATIC(read_pairs_ch)

    trimmed_paired_reshaped_ch = TRIMMOMATIC.out.paired
        .map { sample_id, paired_reads ->
            def forward_paired = paired_reads.find {
                it.name == "${sample_id}_forward_paired.fq.gz"
            }
            def reverse_paired = paired_reads.find {
                it.name == "${sample_id}_reverse_paired.fq.gz"
            }

            if (!forward_paired || !reverse_paired) {
                throw new IllegalStateException(
                    "Could not identify both paired Trimmomatic outputs for ${sample_id}: ${paired_reads}"
                )
            }

            tuple(sample_id, forward_paired, reverse_paired)
        }

    SPADES(TRIMMOMATIC.out.paired)
    REDUCE_CONTIGS(SPADES.out.contigs, params.min_contig_length)

    BACTINSPECTOR(SPADES.out.contigs, batch_name)
    BACTINSPECTOR_BATCH_SUMMARY(
        BACTINSPECTOR.out.tsv.map { sample_id, report -> report }.collect().ifEmpty([]),
        BACTINSPECTOR.out.tsv.map { sample_id, report -> sample_id }.collect().ifEmpty([]),
        expected_sample_ids_ch,
        batch_name
    )
    QUAST(SPADES.out.contigs, batch_name)
    QUAST_BATCH_SUMMARY(
        QUAST.out.tsv.map { sample_id, report -> report }.collect().ifEmpty([]),
        QUAST.out.tsv.map { sample_id, report -> sample_id }.collect().ifEmpty([]),
        expected_sample_ids_ch,
        batch_name
    )

    QUAST_BACTINSPECTOR_FLAGS(
        QUAST_BATCH_SUMMARY.out.tsv,
        QUAST_BATCH_SUMMARY.out.excluded,
        BACTINSPECTOR_BATCH_SUMMARY.out.tsv,
        BACTINSPECTOR_BATCH_SUMMARY.out.excluded,
        expected_sample_ids_ch,
        batch_name,
        Channel.value(file("${projectDir}/bin/quast_bactinspector_flags.py")),
        Channel.value(file("${projectDir}/bin/qc_flagging.py"))
    )

    quast_bactinspector_flagged_sample_ids_ch = QUAST_BACTINSPECTOR_FLAGS.out.flagged_ids
        .splitText()
        .map { it.trim() }
        .filter { it }

    quast_bactinspector_flagged_sample_ids_set_ch =
        quast_bactinspector_flagged_sample_ids_ch
            .toList()
            .map { sample_ids -> sample_ids.toSet() }

    BACTINSPECTOR_FILTERED(REDUCE_CONTIGS.out.clean_contigs, batch_name)
    BACTINSPECTOR_FILTERED_BATCH_SUMMARY(
        BACTINSPECTOR_FILTERED.out.tsv.map { sample_id, report -> report }.collect().ifEmpty([]),
        BACTINSPECTOR_FILTERED.out.tsv.map { sample_id, report -> sample_id }.collect().ifEmpty([]),
        expected_sample_ids_ch,
        batch_name
    )
    QUAST_FILTERED(REDUCE_CONTIGS.out.clean_contigs, batch_name)
    QUAST_FILTERED_BATCH_SUMMARY(
        QUAST_FILTERED.out.tsv.map { sample_id, report -> report }.collect().ifEmpty([]),
        QUAST_FILTERED.out.tsv.map { sample_id, report -> sample_id }.collect().ifEmpty([]),
        expected_sample_ids_ch,
        batch_name
    )

    if (params.kraken2_db) {
        kraken2_db_ch = Channel.value(file(params.kraken2_db, checkIfExists: true))
        KRAKEN2(TRIMMOMATIC.out.paired, kraken2_db_ch)
        KRAKEN2_TAXID_SUMMARY(
            KRAKEN2.out.report,
            Channel.value(file("${projectDir}/bin/summarize_kraken2_report.py")),
            params.target_genus_taxid,
            params.target_species_taxid
        )
        KRAKEN2_BATCH_SUMMARY(
            KRAKEN2.out.report.map { sample_id, report -> report }.collect(),
            KRAKEN2_TAXID_SUMMARY.out.csv.map { sample_id, csv -> csv }.collect(),
            batch_name,
            Channel.value(file("${projectDir}/bin/kraken2_taxa_summary.py")),
            Channel.value(file("${projectDir}/bin/qc_flagging.py"))
        )

        trimmed_paired_reassembly_ch = trimmed_paired_reshaped_ch
            .combine(quast_bactinspector_flagged_sample_ids_set_ch)
            .filter { sample_id, forward_paired, reverse_paired, flagged_sample_ids ->
                flagged_sample_ids.contains(sample_id)
            }
            .map { sample_id, forward_paired, reverse_paired, flagged_sample_ids ->
                tuple(sample_id, forward_paired, reverse_paired)
            }

        kraken2_results_reassembly_ch = KRAKEN2.out.results
            .combine(quast_bactinspector_flagged_sample_ids_set_ch)
            .filter { sample_id, kraken_results, flagged_sample_ids ->
                flagged_sample_ids.contains(sample_id)
            }
            .map { sample_id, kraken_results, flagged_sample_ids ->
                tuple(sample_id, kraken_results)
            }

        NG_READ_LIST(
            kraken2_results_reassembly_ch,
            params.target_genus_taxid,
            params.target_species_taxid
        )

        seqtk_ng_filter_input_ch = trimmed_paired_reassembly_ch
            .join(NG_READ_LIST.out.ng_list)

        SEQTK_NG_FILTER(seqtk_ng_filter_input_ch)

        spades_reassembly_input_ch = SEQTK_NG_FILTER.out.ng_filtered_reads
            .map { sample_id, forward_filtered, reverse_filtered ->
                tuple(sample_id, [forward_filtered, reverse_filtered])
            }

        SPADES_REASSEMBLY(spades_reassembly_input_ch)
        REDUCE_CONTIGS_REASSEMBLY(SPADES_REASSEMBLY.out.contigs, params.min_contig_length)
        reassembly_clean_contigs_ch = REDUCE_CONTIGS_REASSEMBLY.out.clean_contigs
        BACTINSPECTOR_REASSEMBLY(SPADES_REASSEMBLY.out.contigs, batch_name)
        QUAST_REASSEMBLY(SPADES_REASSEMBLY.out.contigs, batch_name)
    } else {
        reassembly_clean_contigs_ch = Channel.empty()
        log.warn 'Kraken2 skipped: provide the required --kraken2_db parameter to run it.'
    }

    assembly_candidates_ch = REDUCE_CONTIGS.out.clean_contigs
        .join(reassembly_clean_contigs_ch, by: 0, remainder: true)

    use_reassembly_for_sample_ids = (params.use_reassembly_for ?: [])
        .collect { it.toString() }
        .toSet()

    assembly_for_downstream_ch = assembly_candidates_ch
        .map { sample_id, original_contigs, reassembled_contigs ->
            if (!original_contigs) {
                throw new IllegalStateException(
                    "Original clean contigs are missing for ${sample_id}"
                )
            }
            if (use_reassembly_for_sample_ids.contains(sample_id)) {
                if (!reassembled_contigs) {
                    throw new IllegalStateException(
                        "Reassembly was requested for ${sample_id}, but no " +
                        "REDUCE_CONTIGS_REASSEMBLY output exists"
                    )
                }
                tuple(sample_id, reassembled_contigs)
            } else {
                tuple(sample_id, original_contigs)
            }
        }

    assembly_selection_rows_ch = assembly_for_downstream_ch
        .map { sample_id, selected_contigs ->
            def promoted = use_reassembly_for_sample_ids.contains(sample_id)
            tuple(sample_id, promoted ? 'reassembled' : 'original', promoted.toString())
        }

    ASSEMBLY_SELECTION(assembly_selection_rows_ch.toList(), batch_name)

    emit:
    trimmed_paired_ch      = TRIMMOMATIC.out.paired
    assembly_for_downstream_ch
    expected_sample_ids_ch
}
