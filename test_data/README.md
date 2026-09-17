# test_data/

Two synthetic paired-end samples, SYNTHETIC and SYNTHETIC_B, simulated
directly from a reference genome (wgsim/dwgsim-style reads). Useful for
smoke-testing the pipeline without real patient data.

## Known limitation: KRAKEN2_BATCH_SUMMARY

Because these reads are simulated with no noise or contamination, Kraken2
classifies every single read - unclassified_pct comes out at exactly 0.00%
for both samples.

bin/kraken2_taxa_summary.py plots Kraken2 QC results on a log scale,
which can't place a point at exactly 0% (log(0) is undefined). So it
raises an error instead of silently misplotting:

  ValueError: Log-scale QC plot requires positive percentages for every sample

Every other stage of the pipeline runs cleanly against these fixtures.
Only this one plotting step fails, and only because of how clean the
simulated reads are - real sequencing data essentially never lands at
exactly 0% unclassified.

This is left as-is for now rather than patched, since loosening the check
would change real behavior in a way that hasn't been validated. Fixing it
properly means regenerating the synthetic reads with a small amount of
realistic noise so unclassified_pct lands above zero.
