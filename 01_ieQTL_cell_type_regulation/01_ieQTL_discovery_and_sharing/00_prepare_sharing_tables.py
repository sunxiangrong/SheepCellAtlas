#!/usr/bin/env python3
"""Prepare auditable tissue-sharing source tables for the ieQTL Extended Data figure.

This script intentionally uses only the Python standard library.  It does not
render figures.  Two sharing definitions are retained separately:

1. tissue-level: a gene-variant signal is shared across tissues, regardless of
   cell type;
2. matched-cell-type: a cell-type-specific gene-variant signal is shared across
   tissues for the same cell type.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def read_gzip_tsv(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, fieldnames, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--merged", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    pair_rows = list(read_tsv(args.summary))
    pair_lookup = {(r["tissue"], r["cell_type"], r["pair"]): r for r in pair_rows}

    tissue_signal_tissues = defaultdict(set)
    tissue_signal_cells = defaultdict(set)
    pair_signals = defaultdict(set)
    pair_matched_signals = defaultdict(set)
    matched_signal_tissues = defaultdict(set)

    for row in read_gzip_tsv(args.merged):
        tissue = row["tissue"]
        cell_type = row["cell_type"]
        pair = row["pair"]
        gene = row["pheno_id"]
        variant = row["variant_id"]
        signal = (gene, variant)
        matched_signal = (cell_type, gene, variant)
        tissue_signal_tissues[signal].add(tissue)
        tissue_signal_cells[signal].add((tissue, cell_type))
        matched_signal_tissues[matched_signal].add(tissue)
        pair_signals[(tissue, cell_type, pair)].add(signal)
        pair_matched_signals[(tissue, cell_type, pair)].add(matched_signal)

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    tissue_combo = defaultdict(lambda: {"n_signals": 0, "n_genes": set()})
    for (gene, variant), tissues in tissue_signal_tissues.items():
        combo = tuple(sorted(tissues))
        tissue_combo[combo]["n_signals"] += 1
        tissue_combo[combo]["n_genes"].add(gene)
    upset_rows = []
    for combo, values in sorted(tissue_combo.items(), key=lambda item: (-item[1]["n_signals"], item[0])):
        upset_rows.append(
            {
                "tissue_set": ";".join(combo),
                "n_tissues": len(combo),
                "n_gene_variant_signals": values["n_signals"],
                "n_genes": len(values["n_genes"]),
            }
        )
    write_tsv(
        out / "tissue_upset_summary.tsv",
        ["tissue_set", "n_tissues", "n_gene_variant_signals", "n_genes"],
        upset_rows,
    )

    all_cells = sorted({(r["tissue"], r["cell_type"]) for r in pair_rows})
    cell_rows = []
    for tissue, cell_type in all_cells:
        pair = next((r for r in pair_rows if r["tissue"] == tissue and r["cell_type"] == cell_type), None)
        if pair is None:
            continue
        pair_key = (tissue, cell_type, pair["pair"])
        tissue_signals = pair_signals.get(pair_key, set())
        matched_signals = pair_matched_signals.get(pair_key, set())
        total = len(tissue_signals)
        shared = sum(len(tissue_signal_tissues[signal]) > 1 for signal in tissue_signals)
        matched_total = len(matched_signals)
        matched_shared = sum(len(matched_signal_tissues[signal]) > 1 for signal in matched_signals)
        cell_rows.append(
            {
                "tissue": tissue,
                "cell_type": cell_type,
                "pair": pair_key[2],
                "n_tissue_level_signals": total,
                "n_tissue_shared_signals": shared,
                "tissue_shared_fraction": f"{shared / total:.10g}" if total else "NA",
                "n_within_tissue_celltype_shared_signals": sum(
                    sum(1 for candidate_tissue, _ in tissue_signal_cells[signal] if candidate_tissue == tissue) > 1
                    for signal in tissue_signals
                ),
                "within_tissue_celltype_shared_fraction": f"{sum(sum(1 for candidate_tissue, _ in tissue_signal_cells[signal] if candidate_tissue == tissue) > 1 for signal in tissue_signals) / total:.10g}" if total else "NA",
                "n_matched_cell_type_signals": matched_total,
                "n_matched_cell_type_shared_signals": matched_shared,
                "matched_cell_type_shared_fraction": f"{matched_shared / matched_total:.10g}" if matched_total else "NA",
            }
        )
    write_tsv(
        out / "tissue_celltype_sharing.tsv",
        [
            "tissue",
            "cell_type",
            "pair",
            "n_tissue_level_signals",
            "n_tissue_shared_signals",
            "tissue_shared_fraction",
            "n_within_tissue_celltype_shared_signals",
            "within_tissue_celltype_shared_fraction",
            "n_matched_cell_type_signals",
            "n_matched_cell_type_shared_signals",
            "matched_cell_type_shared_fraction",
        ],
        cell_rows,
    )

    write_tsv(
        out / "sharing_definition.tsv",
        ["analysis", "signal_unit", "sharing_rule"],
        [
            {
                "analysis": "tissue_upset",
                "signal_unit": "gene-variant",
                "sharing_rule": "same gene-variant observed in more than one tissue, regardless of cell type",
            },
            {
                "analysis": "celltype_heatmap",
                "signal_unit": "cell-type-gene-variant",
                "sharing_rule": "same cell-type-gene-variant observed in more than one tissue",
            },
        ],
    )

    # Gene-level companion tables.  These collapse all significant variants
    # for a gene within each tissue-cell pair to one ieGene, so LD-linked
    # variants do not inflate the displayed sharing counts.
    pair_genes = defaultdict(set)
    tissue_gene_tissues = defaultdict(set)
    tissue_gene_cells = defaultdict(set)
    for row in read_gzip_tsv(args.merged):
        tissue = row["tissue"]
        cell_type = row["cell_type"]
        pair = row["pair"]
        gene = row["pheno_id"]
        pair_genes[(tissue, cell_type, pair)].add(gene)
        tissue_gene_tissues[gene].add(tissue)
        tissue_gene_cells[gene].add((tissue, cell_type))

    gene_combo = defaultdict(set)
    for gene, tissues in tissue_gene_tissues.items():
        gene_combo[tuple(sorted(tissues))].add(gene)
    gene_upset_rows = []
    for combo, genes in sorted(gene_combo.items(), key=lambda item: (-len(item[1]), item[0])):
        gene_upset_rows.append({
            "tissue_set": ";".join(combo),
            "n_tissues": len(combo),
            "n_shared_ieGenes": len(genes),
        })
    write_tsv(
        out / "tissue_upset_ieGene_summary.tsv",
        ["tissue_set", "n_tissues", "n_shared_ieGenes"],
        gene_upset_rows,
    )

    gene_cell_rows = []
    for tissue, cell_type in all_cells:
        pair_row = pair_lookup.get((tissue, cell_type, next(
            r["pair"] for r in pair_rows
            if r["tissue"] == tissue and r["cell_type"] == cell_type
        )))
        pair = pair_row["pair"] if pair_row else ""
        genes = pair_genes.get((tissue, cell_type, pair), set())
        shared_genes = {
            gene for gene in genes
            if sum(1 for candidate_tissue, _ in tissue_gene_cells[gene]
                   if candidate_tissue == tissue) > 1
        }
        gene_cell_rows.append({
            "tissue": tissue,
            "cell_type": cell_type,
            "pair": pair,
            "n_ieGenes": len(genes),
            "n_within_tissue_shared_ieGenes": len(shared_genes),
            "within_tissue_shared_ieGene_fraction": (
                f"{len(shared_genes) / len(genes):.10g}" if genes else "NA"
            ),
        })
    write_tsv(
        out / "tissue_celltype_ieGene_sharing.tsv",
        [
            "tissue",
            "cell_type",
            "pair",
            "n_ieGenes",
            "n_within_tissue_shared_ieGenes",
            "within_tissue_shared_ieGene_fraction",
        ],
        gene_cell_rows,
    )
    write_tsv(
        out / "sharing_definition_ieGene.tsv",
        ["analysis", "signal_unit", "sharing_rule"],
        [
            {
                "analysis": "tissue_upset",
                "signal_unit": "ieGene",
                "sharing_rule": "same ieGene observed in more than one tissue; variants are collapsed within each tissue",
            },
            {
                "analysis": "within_tissue_celltype",
                "signal_unit": "ieGene",
                "sharing_rule": "same ieGene observed in at least two cell types from the same tissue",
            },
        ],
    )

    print(f"pair_rows={len(pair_rows)}")
    print(f"pair_signal_sets={len(pair_signals)}")
    print(f"tissue_signals={len(tissue_signal_tissues)}")
    print(f"matched_celltype_signals={len(matched_signal_tissues)}")
    print(f"upset_combinations={len(upset_rows)}")
    print(f"ieGene_upset_combinations={len(gene_upset_rows)}")


if __name__ == "__main__":
    main()
