# Human disease GWAS mapped to sheep cell contexts

## Intended sequence

1. Obtain each disease's GWAS/COJO outputs and record source, version, build, and coordinate convention.
2. Run `00_prepare_500kb_gwas_loci.py --help` to crop merged per-trait loci and write a processing manifest. This is not raw GWAS harmonization or coordinate conversion.
3. Standardize variant IDs, alleles, and gene mapping; prepare the sheep atlas and disease-specific human reference objects/cell labels.
4. Run the exact scPagwas workflow and summarize scores/FDR by lineage, tissue, and cell type.
5. Run the asthma-specific `04_asthma_lung_airway_case.py` only after the expected summarized association and expression inputs are available. Explicitly select the human ortholog symbol column (`hg19` or `hg38`) to match the reference object.

This folder does not yet contain a portable implementation of every step for all diseases. The asthma script is a case-stage module, not a substitute for GWAS mapping. Coordinate conversion and disease-by-tissue scPagwas job orchestration still need source review and parameterization.
