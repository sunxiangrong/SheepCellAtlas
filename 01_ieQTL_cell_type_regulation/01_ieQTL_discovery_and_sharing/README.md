# ieQTL discovery results and sharing

## Intended sequence

1. Produce bulk-sample cell fractions from the single-cell reference. The study used DWLS; its exact portable command/configuration is not yet in this package.
2. Harmonize genotype, expression, cell fractions, sample covariates, and tissue/cell-type labels; fit genotype × cell-fraction models with the study's OmiGA version. This upstream stage is not included yet.
3. Export a merged significant-ieQTL table and tissue–cell-type pair summary using the thresholds declared by the study.
The analysis-wide stages 1–3 above remain unbundled. The three scripts below are numbered locally from `00` and operate on those upstream outputs:

1. `00_prepare_sharing_tables.py --summary ... --merged ... --out-dir ...` — generate signal-level and ieGene-level sharing tables.
2. `01_plot_detectability.py --help` — plot a compatible precomputed detectability summary.
3. `02_plot_ieGene_sharing.py --merged ... --out-dir ...` — visualize across-tissue and within-tissue cell-type sharing and save the source tables used by the plots.

## Table requirements

The merged ieQTL table is gzip-compressed TSV and needs at least `tissue`, `cell_type`, `pair`, `pheno_id`, and `variant_id`. The pair summary needs `tissue`, `cell_type`, and `pair`. The plot input needs `tissue`, `cell_type`, `pair`, `effective_fraction_sample_size`, and `iegene_detection_proportion` (0–1).

The sharing script is a downstream summarizer; it does not estimate fractions, fit OmiGA, or determine significance thresholds. Record those upstream decisions with the data manifest.
