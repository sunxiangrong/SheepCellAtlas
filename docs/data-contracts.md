# Input and output contracts

This document describes the current portable scripts. It does not imply that missing upstream stages are implemented. Tables use UTF-8 TSV with headers unless the relevant script states otherwise. Controlled or individual-level inputs must stay outside Git.

## 01 · ieQTL genetic regulation

### Sharing table preparation

`01_ieQTL_cell_type_regulation/01_ieQTL_discovery_and_sharing/00_prepare_sharing_tables.py` takes a pair-summary TSV (`tissue`, `cell_type`, `pair`) and merged significant ieQTL gzip TSV (`tissue`, `cell_type`, `pair`, `pheno_id`, `variant_id`). It produces gene-level and variant-level sharing summaries plus the explicit sharing definition.

`01_ieQTL_cell_type_regulation/01_ieQTL_discovery_and_sharing/01_plot_detectability.py` consumes `tissue`, `cell_type`, `pair`, `effective_fraction_sample_size`, and `iegene_detection_proportion` (0–1 scale). These are downstream outputs; DWLS cell fraction estimation and OmiGA model fitting are not in the package.

### Colocalization

`supporting_analyses_not_in_requested_figures/01_GWAS_ieQTL_colocalization/01_run_colocalization.R` expects `chr1.tsv` through `chr26.tsv` inside each trait folder. GWAS columns are `SNP`, `Freq`, `b`, `se`, and `p`; ieQTL columns are `pheno_id`, `variant_id`, `af`, `beta_g2`, `beta_se_g2`, and `pval_g2`.

`02_plot_candidate_loci.py` expects a candidate manifest with `trait`, `tissue`, `cell_type`, `gene`, `lead_snp`, `pph3`, and `pph4`, plus regional GWAS and ieQTL files. `03_plot_ear_length_case.py` additionally needs a PLINK genotype prefix and executable to calculate LD. Confirm assembly, chromosome labels, alleles, and variant IDs before use.

Bulk eQTL–ieQTL comparison is a distinct branch and is not currently represented by a full portable script sequence. Do not treat the GWAS–ieQTL outputs as a bulk-QTL comparison.

## 02 · Complex-trait cell mapping

The result-summary plots expect a result root arranged as `<tissue>/<trait>/Merged_celltype_summary.csv`, plus a lineage-map TSV with `cell_type` and `cell_lineage`. Summary tables require `celltype`, `mean_TRS`, `merged_BG_adjp`, and `n_cells`.

The context-ranking script expects context-specific `pcc_top1000_<context>.tsv` tables with `gene` and `PCC_clean`. Input trait list, phenotype units, reference assembly, scPagwas version, thresholds, and mapping configuration must be supplied with the upstream data manifest; those upstream stages are not currently bundled.

## 03 · Cross-species disease analysis

The expression-bridge scripts and asthma case require one-to-one ortholog records with `Human_homology_type`, `Sheep_gene_name`, `Human_gene_name_hg19`, and optionally `Human_gene_name_hg38`. Explicitly choose the matching gene-symbol field with `--human-ortholog-column`; the scripts should not silently guess the build. The Seurat exporter writes cell-level group means, group cell counts, and session information for an explicitly selected group column and RNA layer; this is not sample-aware raw-count pseudobulk. Sheep expression tables begin with `gene`; human `.h5ad` inputs need the requested grouping column and declared expression layer. Cell labels must already be harmonized.

The human GWAS locus-window preparation script expects already processed and merged per-trait GWAS/COJO loci. It writes locus windows and a processing manifest; it does not harmonize raw summary statistics, convert genome builds, or align alleles. Confirm build and coordinate convention before use.

The psoriasis summary expects `candidate_pcc_with_stable_ortholog_flag.tsv.gz` and `candidate_stable_fraction_enrichment_summary.tsv`. The dermatitis GO plotting script expects precomputed conserved-DEG and GO-enrichment tables. The asthma case requires curated human disease-association results, a sheep lung pseudobulk, a human airway object, and the ortholog table.

For each disease, record source/version, assembly, coordinate conversion, expression layer, pseudobulk grouping, ortholog release, cell-label crosswalk, enrichment background, and thresholds. Disease GWAS-to-cell mapping and enrichment generation are not fully packaged yet.
