# SheepCA genetics analysis code

The repository follows three manuscript analysis topics. Directories are organized by analysis module and processing order, not by figure number. This approach is informed by the thematic, numbered workflows in the [FarmGTEx Cattle Cell Atlas pipeline](https://github.com/FarmGTEx/CattleCellAtlas_pipeline_V0).

```text
01_ieQTL_cell_type_regulation/
02_complex_trait_cell_mapping/
03_cross_species_disease_analysis/
supporting_analyses_not_in_requested_figures/
```

The Chinese README in each module is the current scope and script map. Missing upstream procedures are labeled as gaps; plotting scripts are not presented as substitutes for statistical analyses. ieQTL refers to genotype-by-inferred-cell-fraction interaction in bulk samples, not direct single-cell eQTL.

This is a code-release candidate, not a complete raw-data-to-result reproduction. It contains no restricted genotypes, expression matrices, single-cell objects, or individual-level metadata. See `PUBLICATION_MANIFEST.tsv` for script status and known limitations, and `docs/release-checklist.md` before publication.
