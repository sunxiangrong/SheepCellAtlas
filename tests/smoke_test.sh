#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
python -m py_compile "$repo_root"/01_ieQTL_cell_type_regulation/01_ieQTL_discovery_and_sharing/*.py "$repo_root"/supporting_analyses_not_in_requested_figures/01_GWAS_ieQTL_colocalization/*.py "$repo_root"/supporting_analyses_not_in_requested_figures/02_regulatory_annotation/*.py "$repo_root"/02_complex_trait_cell_mapping/03_candidate_genes_and_pathways/*.py 02_complex_trait_cell_mapping/01_global_trait_association/*.py "$repo_root"/02_complex_trait_cell_mapping/04_tissue_level_extensions/*.py "$repo_root"/03_cross_species_disease_analysis/01_human_GWAS_to_sheep_cells/*.py "$repo_root"/03_cross_species_disease_analysis/02_conserved_expression_and_enrichment/*.py

python "$repo_root"/01_ieQTL_cell_type_regulation/01_ieQTL_discovery_and_sharing/01_plot_detectability.py --help >/dev/null
python "$repo_root"/01_ieQTL_cell_type_regulation/01_ieQTL_discovery_and_sharing/02_plot_ieGene_sharing.py --help >/dev/null
python "$repo_root"/supporting_analyses_not_in_requested_figures/02_regulatory_annotation/00_e5_matched_background_analysis.py --help >/dev/null
python "$repo_root"/supporting_analyses_not_in_requested_figures/01_GWAS_ieQTL_colocalization/02_plot_candidate_loci.py --help >/dev/null
python "$repo_root"/supporting_analyses_not_in_requested_figures/01_GWAS_ieQTL_colocalization/03_plot_ear_length_case.py --help >/dev/null
python "$repo_root"/02_complex_trait_cell_mapping/04_tissue_level_extensions/00_plot_trait_landscape.py --help >/dev/null
python "$repo_root"/02_complex_trait_cell_mapping/01_global_trait_association/01_plot_lineage_associations.py --help >/dev/null
python "$repo_root"/02_complex_trait_cell_mapping/03_candidate_genes_and_pathways/00_rank_candidate_contexts.py --help >/dev/null
python "$repo_root"/03_cross_species_disease_analysis/01_human_GWAS_to_sheep_cells/00_prepare_500kb_gwas_loci.py --help >/dev/null
python "$repo_root"/03_cross_species_disease_analysis/02_conserved_expression_and_enrichment/02_tissue_expression_bridge.py --help >/dev/null
python "$repo_root"/03_cross_species_disease_analysis/02_conserved_expression_and_enrichment/03_celltype_expression_bridge.py --help >/dev/null
python "$repo_root"/03_cross_species_disease_analysis/02_conserved_expression_and_enrichment/06_plot_psoriasis_orthologs.py --help >/dev/null
python "$repo_root"/03_cross_species_disease_analysis/01_human_GWAS_to_sheep_cells/04_asthma_lung_airway_case.py --help >/dev/null
