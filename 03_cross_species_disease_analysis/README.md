# 03 · Cross-species analysis of human disease

Workflow: map human disease GWAS to sheep cell contexts → inspect associated expression programs → compare baseline expression of human/sheep disease genes → interpret conserved genes with functional enrichment.

Keep three evidence streams distinct: TRS/cell association; PCC enrichment of stable-expression orthologs; and conserved disease DEGs with GO analysis. Current code includes reusable expression bridges and partial case-result plotting, but not the full disease mapping and figure-specific DEG/GO generation chain. Atopic dermatitis is used for expression and functional analysis only; no AD–TRS result is claimed. The specified functional enrichment is GO, not KEGG.
