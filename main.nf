nextflow.enable.dsl=2

include { QC_ASSEMBLY } from './workflows/qc_assembly.nf'
include { TYPING_AMR } from './workflows/typing_amr.nf'
include { PHYLOGENY } from './workflows/phylogeny.nf'

workflow {
    // TODO: wire channels once modules are ported
}
