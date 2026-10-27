# GWAS–ieQTL colocalization

## Intended sequence

1. Harmonize GWAS and ieQTL variant IDs, alleles, genome build, allele frequencies, and cis windows.
2. Run `01_run_colocalization.R` for per-trait/cell-context colocalization.
3. Prepare a candidate manifest and run `02_plot_candidate_loci.py` for selected regions.
4. Use `03_plot_ear_length_case.py` for the ear-length `STXBP1` example, including LD from the supplied PLINK genotype files.

The R script expects chromosome-partitioned GWAS input (`chr1.tsv`–`chr26.tsv`) and the columns listed in [`../../docs/data-contracts.md`](../../docs/data-contracts.md). The colocalization region is determined by the upstream input; this script does not itself enforce a ±500-kb window. The minimum shared-SNP count is supplied as an argument. Use a fresh output directory; failed colocalizations are recorded in `coloc_errors.tsv`. `00_compute_locus_ld.py` exports PLINK 1.9 LD separately, with the command and version recorded. The plot scripts are result visualization, not a substitute for data harmonization or the colocalization calculation.
