# N. gonorrhoeae WGS pipeline — manual scripts

16 standalone scripts, one per pipeline step. Each script has an "EDIT THESE"
block at the top — set the paths and sample name there, then run the script.
No shared config file, no dependencies between scripts other than each one
reading the output of an earlier step.

> **Scope note**: these scripts mirror the core WGS pipeline (QC through
> phylogeny) as standalone bash, useful for understanding or debugging any
> pipeline step without Nextflow. They do not yet cover the reassembly/
> decontamination branch or SensiTyper — those are Nextflow-only for now.

---

## Order

| # | Script | Input needed |
|---|---|---|
| 1 | `01_qc.sh` | raw FASTQs |
| 2 | `02_trimmomatic.sh` | raw FASTQs |
| 3 | `03_kraken2.sh` | output of step 2 |
| 4 | `04_spades.sh` | output of step 2 |
| 5 | `05_bactinspector.sh` | output of step 4 |
| 6 | `06_quast.sh` | output of step 4, all samples at once |
| 7 | `07_ariba.sh` | output of step 2 |
| 8 | `08_pyngost.sh` | output of step 4, all samples at once |
| 9 | `09_snippy_reference.sh` | output of step 2, run per sample |
| 10 | `10_snippy_core.sh` | output of step 9, once per batch |
| 11 | `11_snippy_assembly.sh` | output of steps 2 and 4, run per sample |
| 12 | `12_het_sites_filter.sh` | output of step 9 or 11, run per sample |
| 13 | `13_het_sites_panel.sh` | output of step 12, once per batch |
| 14 | `14_het_sites_requery.sh` | output of step 13, run per sample |
| 15 | `15_prokka.sh` | output of step 4 |
| 16 | `16_gubbins_iqtree.sh` | output of step 10 |

Steps 9 and 11 are two different analyses (different reference, different
thresholds) — run both if you want the full het-site QC.

Steps 12–14 run once for whichever branch you're QC'ing (reference-mapped or
assembly-mapped) — point `RAW_VCF` at the right one.

## Before running

- Each tool needs its own environment active (conda env or module) — activate
  it before running that script.
- Steps 6, 8, 10, and 13 are batch-wide: run them once, after the per-sample
  steps are done for every sample.
