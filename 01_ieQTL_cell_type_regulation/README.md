# 01 · Cell-type interaction eQTLs

This module covers genotype-by-inferred-cell-fraction interaction in bulk samples, followed by comparison with bulk eQTLs. It does not describe direct single-cell eQTLs.

The reusable scripts in `01_ieQTL_discovery_and_sharing/` prepare sharing summaries and plot precomputed results. The `02_bulk_eQTL_ieQTL_comparison/` folder records the comparison scope; its complete portable workflow is still missing. DWLS fraction estimation, the study's OmiGA model run, and the selected nonduplicative case scripts also remain to be added.

The current audited counts are 980 tissue–cell-type–gene ieQTL records, 879 unique ieGenes, and 336 records whose corresponding-tissue bulk eGene did not meet `qval_g1 < 0.05`. This means below the bulk eGene significance threshold, not absence of a gene or bulk signal. GWAS–ieQTL colocalization and E5 annotation are kept under `supporting_analyses_not_in_requested_figures/`.
