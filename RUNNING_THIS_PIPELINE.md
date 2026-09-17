## Requirements

- Nextflow (tested on 26.04.6)
- Apptainer or Singularity, or Docker (Apptainer is the primary target;
  Docker support exists but isn't the main path)
- Enough memory for SPAdes — 64 GB by default, set with --spades_mem

## Input data layout

Paired-end FASTQ files named `<sample>_R1_001.fastq.gz` and
`<sample>_R2_001.fastq.gz` in one directory, passed with --input
(defaults to test_data/).

## Running locally

Minimum command, using the bundled synthetic test data and no external
databases:

    nextflow run main.nf -profile singularity

With no databases set, Kraken2, ARIBA, pyngoST, and the Snippy branches
are skipped with a warning rather than failing. To run everything:

    nextflow run main.nf -profile singularity \
        --kraken2_db /path/to/kraken2_db \
        --ariba_db /path/to/ariba_db.tar.gz \
        --pyngost_db /path/to/pyngost_allelesDB \
        --snippy_reference /path/to/reference.fasta \
        --snippy_mask /path/to/mask.bed

The reference FASTA needs a BWA index alongside it (.amb .ann .bwt .pac
.sa). --run_phylogeny turns on the Gubbins recombination-masking step,
once Snippy-core has run; --run_assembly_mapped_snippy turns on the
assembly-mapped Snippy branch (its heterozygous-site QC isn't built yet,
so only the variant calling itself runs).

## Running with containers

Every step runs in a pinned container. Most are pulled automatically
from public registries the first time they're needed. Two are built
locally instead, since no ready-made image exists for them:

    apptainer build containers/bactinspector/bactinspector.sif containers/bactinspector/bactinspector.def
    apptainer build containers/sensitype/sensitype.sif containers/sensitype/sensitype.def

Both .sif files are build artifacts, not source — they're gitignored
and need rebuilding on a new machine.

## Interpreting output

Results are written under results/, organized by pipeline stage and by
batch name (defaults to the Nextflow run name if you don't set
--batch_name explicitly — pass one if you want stable, comparable output
directories across runs). Each per-batch stage also writes an
"excluded/failed samples" file, listing any expected sample that didn't
produce output, with a reason.

## Running steps manually, without Nextflow

See docs/manual_reference/ — one script per step, for debugging outside
the pipeline.
