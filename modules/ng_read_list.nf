/*
 * Consume Kraken2's raw per-read classification output (not its summary
 * report) and produce the read-ID list used for the targeted-reassembly
 * filter. Target taxids are parameterized (default: NCBI taxid 482 =
 * genus Neisseria, 485 = species N. gonorrhoeae) so a different organism
 * can supply its own without editing this process.
 */
process NG_READ_LIST {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'debian:12.11-slim'
    publishDir 'results/ng_read_list', mode: 'copy'

    input:
    tuple val(sample_id), path(kraken_results)
    val target_genus_taxid
    val target_species_taxid

    output:
    tuple val(sample_id), path("${sample_id}_NG.lst"), emit: ng_list

    script:
    """
    awk '\$1=="C" && (\$3==${target_genus_taxid} || \$3==${target_species_taxid}) {print \$2}' ${kraken_results} > ${sample_id}_NG.lst
    """
}
