# Matched-background regulatory annotation

`00_e5_matched_background_analysis.py` compares tissue-specific FAANG E5 enhancer overlap among PPH3-high and PPH4-high ieQTL leads against a matched cis-testing background.

## Run order

1. Prepare the lead table with `group` (`PPH3-high` or `PPH4-high`), `tissue`, `cell_type`, `phenotype2`, and `variant_id` columns.
2. Supply the per-tissue/cell-type nominal ieQTL files and the matching-assembly E5 BED.
3. Run `python 00_e5_matched_background_analysis.py --help`, then provide a new `--out-dir`; set and retain the random seed and permutation count.
4. Archive the annotated lead table, sampled background distributions, summary table, and plots together.

The control pool is drawn from the same gene's nominal ieQTL variants and matched on allele frequency and absolute cis distance according to the score implemented in the script. This is an enrichment comparison under that matching design, not a genome-wide causal test. Use the assembly and BED interval convention that match the variants.
