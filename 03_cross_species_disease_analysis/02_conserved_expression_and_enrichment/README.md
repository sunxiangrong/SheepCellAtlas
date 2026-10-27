# Conserved expression and enrichment

## Intended sequence

1. Run `00_export_seurat_group_mean.R <seurat_rds> <group_column> <layer> <out_tsv>` to export mean expression, group cell counts, and session information from an annotated Seurat object.
2. Prepare matched human tissue and cell-type expression inputs, documenting the expression layer and aggregation rule.
3. Restrict both species to the selected one-to-one ortholog column and harmonize gene IDs.
4. Run `02_tissue_expression_bridge.py` and `03_celltype_expression_bridge.py` to produce matched cross-species expression summaries. Select `--human-ortholog-column` explicitly to match the human gene symbols in the H5AD.
5. Define conserved disease-associated genes and run the declared enrichment test with its tested-gene background (upstream generation not bundled).
6. Use `05_plot_ad_go_enrichment.R` and `06_plot_psoriasis_orthologs.py` to visualize supplied result tables.

The export calculates a cell-level group mean; it is not raw-count, sample-aware pseudobulk and does not provide biological-replicate estimates. The bridge scripts do not choose orthologs, map cell labels, or perform complete disease DEG/enrichment analysis automatically. Record the ortholog release, aggregation rule, expression layer, and stable-expression criterion.
