/*
 * NG-specific by design: SensiTyper is a tool built specifically for
 * N. gonorrhoeae AMR resistance-profile typing (penA mosaic-allele calling,
 * NG-specific antibiotic panel below). Adapting this pipeline for another
 * organism means replacing this whole module with that organism's own
 * resistance-typing tool, not adjusting its parameters.
 *
 * SensiTyper's sensitype analysis is batch-wide: the ARIBA summary CSV names
 * report_complete.tsv files below the accompanying reports/ directory.
 * Keep both inputs staged together so those relative paths remain resolvable.
 *
 * This process deliberately uses Nextflow's default terminating error strategy.
 * Unlike per-sample processes, there is no independent sample result to retain
 * if this single batch task fails; ignoring it would silently lose the complete
 * resistance-profile report.
 */
process SENSITYPE {
    tag "batch"
    container "${projectDir}/containers/sensitype/sensitype.sif"
    publishDir 'results/sensitype', mode: 'copy'

    input:
    path ariba_summary_csv
    path ariba_reports

    output:
    path "sensiscript_results.tsv", emit: tsv
    path "sensiscript_results.html", emit: html

    script:
    """
    python3 /opt/sensityper/sensiscript_v2.6.py \
        --input_AMRtable ${ariba_summary_csv} \
        --outfile sensiscript_results.tsv \
        --database /opt/sensityper/sensitype.db \
        --pena /opt/sensityper/sensitype.penA.db \
        --antibiotics ceftriaxone,azithromycin,ciprofloxacin,spectinomycin,penicillin,tetracycline,zoliflodacin

    test -s sensiscript_results.tsv
    test -s sensiscript_results.html
    """
}
