# ng-surveillance

A Nextflow pipeline for whole-genome sequencing analysis of a bacterial
organism: quality control, assembly, typing, AMR detection, and
phylogenetics.

It was built for Neisseria gonorrhoeae — genus/species taxid filtering,
the AMR/typing databases, and the reference genome are all organism-
specific — but the QC, assembly, and phylogenetics stages are generic,
and it can be adapted to other organisms by swapping those inputs and
databases.

## What it does

- Read QC and trimming (FastQC, MultiQC, Trimmomatic)
- Taxonomic screening (Kraken2, BactInspector) and assembly QC (SPAdes,
  QUAST), with automatic flagging of samples that need reassembly
- Genome annotation (Prokka)
- AMR gene detection and sequence typing (ARIBA, pyngoST, SensiTyper)
- Reference-mapped variant calling with heterozygous-site QC (Snippy),
  and an assembly-mapped branch for samples without a shared reference
- Recombination-masked core-genome alignment (Snippy-core, Gubbins),
  behind an opt-in flag

## How this compares to IMMense

See docs/IMMENSE_COMPARISON.md. Short version: this pipeline follows the
same engineering standards (Nextflow, pinned containers, an audit trail
of what produced what) but swaps in organism-specific tools where a
general-purpose pipeline doesn't cover what this project needs.

## Quickstart

Requires Nextflow and Apptainer (or Singularity). Two containers are
built locally from definition files before first use:

    apptainer build containers/bactinspector/bactinspector.sif containers/bactinspector/bactinspector.def
    apptainer build containers/sensitype/sensitype.sif containers/sensitype/sensitype.def

Then, against the bundled synthetic test data:

    nextflow run main.nf -profile singularity

This runs QC, assembly, and annotation with no extra setup. AMR typing
and phylogenetics need real databases — see RUNNING_THIS_PIPELINE.md.

## Parameters

See RUNNING_THIS_PIPELINE.md for the full list, and nextflow_schema.json
for the machine-readable version.

## Manual reference

docs/manual_reference/ has one standalone script per pipeline step, for
running or debugging a step outside Nextflow. It doesn't cover the
reassembly branch or SensiTyper.

## Citing

See CITATION.cff.
