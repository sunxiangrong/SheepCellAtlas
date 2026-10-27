# GitHub release plan

## Release boundary

Publish this repository as the **analysis-code companion** to the manuscript. It should contain scripts, input schemas, environment definition, provenance, and lightweight tests. It must not contain raw genotype or expression data, human single-cell objects, individual/sample metadata, unpublished result tables, or rendered manuscript figures.

## Recommended sequence

1. **Scientific freeze**: resolve the analysis gaps in `release-checklist.md` — upstream ieQTL modeling, bulk-eQTL/ieQTL comparison, complex-trait and disease mapping configuration, and virtual-knockout provenance.
2. **Data release plan**: assign every input in `PUBLICATION_MANIFEST.tsv` to either an accession/DOI, a controlled-access request route, or an explicit non-release reason. Add checksums and versions for the released derived tables.
3. **Metadata**: replace placeholders, add `LICENSE`, create `CITATION.cff`, state the manuscript version, and add the final data and article DOI when available.
4. **Private GitHub review**: create a private repository from this folder, run the CI workflow, inspect the complete file list, and have coauthors approve the data-access statement and license.
5. **Public release**: change visibility only after the safety and scientific checks pass; create a semantic release tag such as `v1.0.0` matching the manuscript submission/revision state.
6. **Archive**: connect the tagged GitHub release to Zenodo or the journal-approved archive and cite the archived DOI in the paper.

## Suggested repository identity

- Repository name: `sheepca-genetic-regulation`
- Short description: `Reproducible code for Sheep Cell Atlas ieQTL, complex-trait, and cross-species disease analyses.`
- Initial release: `v1.0.0` only after the workflow and data boundaries above are frozen.

## Publication-host rule

Prepare and review on the Mac, but use the established approved publishing route for the eventual GitHub write. Do not push controlled data or credentials from any host.
