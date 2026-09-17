nextflow.enable.dsl = 2

include { QC_ASSEMBLY } from './workflows/qc_assembly'
include { TYPING_AMR } from './workflows/typing_amr'
include { PHYLOGENY } from './workflows/phylogeny'

/*
 * ng-surveillance — N. gonorrhoeae WGS AMR surveillance pipeline.
 * Parameter defaults live in nextflow.config, not here — see its params {}
 * block. A proper nextflow_schema.json pass (documenting each param,
 * marking required vs. optional) is still pending.
 */

workflow {
    batch_name = params.batch_name ?: workflow.runName

    read_pairs_ch = Channel
        .fromFilePairs(
            "${params.input}/*_R{1,2}_001.fastq.gz",
            checkIfExists: true,
            flat: false
        )

    QC_ASSEMBLY(read_pairs_ch, batch_name)

    TYPING_AMR(
        QC_ASSEMBLY.out.assembly_for_downstream_ch,
        QC_ASSEMBLY.out.trimmed_paired_ch,
        batch_name
    )

    PHYLOGENY(
        QC_ASSEMBLY.out.trimmed_paired_ch,
        QC_ASSEMBLY.out.assembly_for_downstream_ch,
        QC_ASSEMBLY.out.expected_sample_ids_ch,
        batch_name
    )
}
