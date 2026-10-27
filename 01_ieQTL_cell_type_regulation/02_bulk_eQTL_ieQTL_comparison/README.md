# Bulk eQTL–ieQTL comparison

This branch compares molecular-QTL evidence in bulk tissue with cell-type interaction results. It is **not** the GWAS–ieQTL colocalization branch.

## Separate outputs to retain

1. Gene-level overlap/sharing table, with the precise unit (gene × tissue or gene × tissue–cell-type pair) stated.
2. Matched locus-level bulk-QTL and ieQTL association tables for any colocalization or regional plot.
3. Effect-size correlation table, with the compared coefficients (for example, bulk `g1` vs ieQTL `g1` or `g2`) explicit.
4. Regulatory-annotation table, with assembly, annotation source/release, interval convention, and overlap rule recorded.

The original projects contain example scripts and case assets, but the current public bundle does not yet contain one portable, validated sequence that builds all four outputs from raw/source tables. The exact study-selected correlation, overlap thresholds, and complete case inputs must be confirmed before this branch can be represented as fully reproducible code.

The input contract is defined in [`../../docs/data-contracts.md`](../../docs/data-contracts.md). Do not label a matched bulk-QTL track as a three-way GWAS–bulk-QTL–ieQTL colocalization.
