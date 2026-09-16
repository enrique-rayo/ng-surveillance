process SPADES {
    tag "$sample_id"
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/spades:4.3.0--hde4eca7_0'
    publishDir {
        task.process.contains('REASSEMBLY') ? 'results/spades_reassembly' : 'results'
    }, mode: 'copy', saveAs: { filename ->
        task.process.contains('REASSEMBLY')
            ? filename.replaceFirst('^spades/', '')
            : filename
    }
    // A 121 GiB host cannot safely run two assemblies at the 64 GB default cap.
    maxForks 1

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("spades/${sample_id}/contigs.fasta"), emit: contigs
    tuple val(sample_id), path("spades/${sample_id}"), emit: assembly

    script:
    """
    spades.py \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        --careful \
        -o spades/${sample_id} \
        -t ${task.cpus} \
        -m ${task.memory.toGiga()}
    """
}
