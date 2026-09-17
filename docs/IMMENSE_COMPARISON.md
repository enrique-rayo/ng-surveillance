# How this compares to IMMense

IMMense (https://gitlab.uzh.ch/appliedmicrobiologyresearch/immense) is the
ISO-accredited pipeline used by the Swiss Pathogen Surveillance Platform
(SPSP). It's a general-purpose bacterial genomics pipeline. This one is
built for a single organism, so several tools differ.

| Step | IMMense | This pipeline | Why |
|---|---|---|---|
| Read QC | FastQC, MultiQC | FastQC, MultiQC | Same |
| Trimming | Trimmomatic | Trimmomatic | Same |
| Assembly | Unicycler | SPAdes | Both are short-read de novo assemblers; SPAdes performs well on small, GC-uniform genomes |
| Assembly QC | QUAST | QUAST | Same |
| Completeness | BUSCO, CheckM | not used | Not needed for a small, well-characterised genome |
| Contamination / taxonomy | GTDB-tk, MetaPhlAn4, 16S/rMLST BLAST | Kraken2, BactInspector | This pipeline targets one known organism, so it needs sensitive detection of related species rather than open-ended classification |
| Annotation | Bakta | Prokka | Functionally similar; Prokka fits the downstream typing tools better |
| AMR / typing | abricate, AMRFinderPlus | ARIBA, pyngoST, SensiTyper | IMMense's generic AMR screening doesn't cover the organism-specific typing schemes (MLST, NG-STAR, mosaic gene calling) this pipeline needs |
| Phylogenetics | not in scope | Snippy, Gubbins | Added for outbreak/cluster detection |

## What matching IMMense's standard actually means here

Not reusing every tool — several of the substitutions above are deliberate.
It means:

- Nextflow, same as IMMense
- every step runs in a pinned container
- every run leaves a record of which tool version and parameters produced
  which output
- per-batch MultiQC reporting

A profile for running on SPSP's own infrastructure is left for later — it
needs testing against real SPSP compute, which hasn't happened yet.
