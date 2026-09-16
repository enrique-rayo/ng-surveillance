/*
 * NG-specific by design. Unlike qc_assembly.nf, this sub-workflow's tools
 * (SensiTyper, pyngoST/NG-STAR) are built specifically for N. gonorrhoeae —
 * there is no generic parameterization that makes AMR/resistance typing
 * organism-agnostic. Adapting this pipeline for another organism means
 * replacing this sub-workflow's modules entirely, not tuning their params.
 * See per-module comments for specifics.
 */
include { PROKKA } from '../modules/prokka'
include { ARIBA_RUN; ARIBA_SUMMARY } from '../modules/ariba'
include { SENSITYPE } from '../modules/sensitype'
include { PYNGOST } from '../modules/pyngost'

workflow TYPING_AMR {
    take:
    assembly_for_downstream_ch   // tuple(sample_id, clean_contigs)
    trimmed_paired_ch            // tuple(sample_id, [forward_paired, reverse_paired])
    batch_name                   // val

    main:
    qc_control_pattern = java.util.regex.Pattern.compile(
        params.qc_control_sample_pattern ?: '(?i)^NGIVPZEC\\d+(?:_S\\d+)?$'
    )

    PROKKA(assembly_for_downstream_ch, params.prokka_genus ?: '')

    if (params.ariba_db) {
        ariba_db_ch = Channel.value(file(params.ariba_db, checkIfExists: true))
        ARIBA_RUN(trimmed_paired_ch, ariba_db_ch)

        ariba_summary_results_ch = ARIBA_RUN.out.results
            .filter { sample_id, result_dir ->
                !qc_control_pattern.matcher(sample_id).matches()
            }

        ARIBA_SUMMARY(
            ariba_summary_results_ch
                .map { sample_id, result_dir -> result_dir }
                .collect()
        )
        SENSITYPE(ARIBA_SUMMARY.out.csv, ARIBA_SUMMARY.out.reports)
    } else {
        log.warn 'ARIBA skipped: provide the required --ariba_db parameter (prepared-database .tar.gz archive) to run it.'
    }

    if (params.pyngost_db) {
        pyngost_db_ch = Channel.value(file(params.pyngost_db, checkIfExists: true))
        pyngost_exclude_sample_ids = (params.pyngost_exclude_samples ?: '')
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it }
            .toSet()

        pyngost_contigs_ch = assembly_for_downstream_ch
            .filter { sample_id, contigs -> !pyngost_exclude_sample_ids.contains(sample_id) }

        pyngost_excluded_samples_ch = assembly_for_downstream_ch
            .filter { sample_id, contigs -> pyngost_exclude_sample_ids.contains(sample_id) }
            .map { sample_id, contigs -> sample_id }
            .toList()

        PYNGOST(
            pyngost_contigs_ch.map { sample_id, contigs -> contigs }.collect(),
            pyngost_db_ch,
            batch_name,
            pyngost_excluded_samples_ch
        )
    } else {
        log.warn 'pyngoST skipped: provide the required --pyngost_db parameter (allelesDB directory) to run it.'
    }
}
