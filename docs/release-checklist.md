# Public release checklist

## Scientific closure

- [ ] Recover and document the exact bulk-eQTL versus ieQTL comparison definition and source tables; do not conflate this with GWAS–ieQTL colocalization.
- [ ] Recover the virtual-knockout inputs, target selection, software version, and settings, or keep it explicitly marked as unavailable.
- [ ] Recover the authoritative upstream psoriasis/dermatitis disease-GWAS-to-cell mapping scripts and verify the final cross-species result set.
- [ ] Add the missing cell-fraction estimation and ieQTL model-fitting steps if those are intended to be part of the public end-to-end workflow.

## Safety and provenance

- [ ] Do not commit raw genotype, expression, single-cell, sample, or individual-level metadata.
- [ ] Scan staged files for absolute paths, credentials, private inputs, and unpublished result assets.
- [ ] Ship a versioned data manifest with checksums, package versions, expression layer, ortholog release, and cell-label mapping.
- [ ] Replace all data-access, author, citation, and DOI placeholders.
- [ ] Select a software license. None is supplied here because it is a legal authorship decision.

## GitHub sequence

1. Create the repository from this folder only after the checklist passes.
2. Add `LICENSE`, `CITATION.cff`, a release tag, and manuscript/data DOI once available.
3. Push through the established approved GitHub publication route; do not upload raw data.
4. Archive the tagged release to Zenodo or the journal-approved repository and add its DOI to the paper.
