#!/usr/bin/env python3
"""Create ieGene-level across-tissue and within-tissue UpSet figures."""

from __future__ import annotations

import argparse
import gzip
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.0,
    "axes.linewidth": 0.65,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

TISSUE_COLORS = {
    "Adipose": "#fdedd4", "Cecum": "#d4a66b", "Cerebellum": "#e5d023",
    "Cerebral cortex": "#f1eb3a", "Cerebrum": "#ecd9a0", "Colon": "#eabb63",
    "Duodenum": "#e39962", "Heart": "#e8401f", "Hypothalamus": "#7e5d26",
    "Ileum": "#c8923b", "Liver": "#aa8a89", "Longissimus muscle": "#fadbae",
    "Muscle": "#fadbae", "Oviduct": "#c090a6", "Perirenal adipose": "#fdedd4",
    "Pituitary": "#e5ce5c", "Skin": "#e3e5e5", "Spleen": "#902731",
}


def display(value: str) -> str:
    return value.replace("_", " ").replace("plus", "+")


def load_membership(path: Path):
    tissue_gene_tissues = defaultdict(set)
    tissue_gene_cells = defaultdict(lambda: defaultdict(set))
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = pd.read_csv(handle, sep="\t", usecols=["tissue", "cell_type", "pheno_id"])
        for row in reader.itertuples(index=False):
            tissue = str(row.tissue)
            cell = str(row.cell_type)
            gene = str(row.pheno_id)
            tissue_gene_tissues[gene].add(tissue)
            tissue_gene_cells[tissue][gene].add(cell)
    return tissue_gene_tissues, tissue_gene_cells


def build_across(tissue_gene_tissues):
    combos = defaultdict(set)
    for gene, tissues in tissue_gene_tissues.items():
        combos[tuple(sorted(tissues))].add(gene)
    rows = []
    for combo, genes in sorted(combos.items(), key=lambda x: (-len(x[1]), x[0])):
        rows.append({
            "tissue_set": ";".join(combo),
            "n_tissues": len(combo),
            "n_shared_ieGenes": len(genes),
        })
    return pd.DataFrame(rows)


def build_within(tissue_gene_cells):
    rows = []
    for tissue, gene_cells in sorted(tissue_gene_cells.items()):
        combos = defaultdict(set)
        for gene, cells in gene_cells.items():
            if len(cells) >= 2:
                combos[tuple(sorted(cells))].add(gene)
        for combo, genes in sorted(combos.items(), key=lambda x: (-len(x[1]), x[0])):
            rows.append({
                "tissue": tissue,
                "cell_type_set": ";".join(combo),
                "n_cell_types": len(combo),
                "n_shared_ieGenes": len(genes),
            })
    return pd.DataFrame(rows)


def draw_upset(ax_bar, ax_matrix, rows, set_col, title, subtitle, set_label, color="#496b87"):
    rows = rows[rows["n_shared_ieGenes"] > 0].copy()
    rows = rows.sort_values("n_shared_ieGenes", ascending=False).reset_index(drop=True)
    sets = sorted({item for value in rows[set_col] for item in value.split(";")})
    matrix = np.zeros((len(rows), len(sets)), dtype=int)
    for i, value in enumerate(rows[set_col]):
        for item in value.split(";"):
            matrix[i, sets.index(item)] = 1

    x = np.arange(len(rows))
    ax_bar.bar(x, rows["n_shared_ieGenes"], color=color, width=0.78, edgecolor="none")
    ax_bar.set_yscale("log")
    ax_bar.set_ylabel("Shared ieGenes\n(log scale)", labelpad=3)
    ax_bar.set_xlim(-0.6, len(rows) - 0.4)
    ax_bar.set_xticks([])
    ax_bar.tick_params(length=2, pad=2)
    ax_bar.set_title(title, loc="left", fontsize=8.5, fontweight="bold", pad=12)
    ax_bar.text(0, 1.02, subtitle, transform=ax_bar.transAxes, fontsize=6.0, color="#52606d", va="bottom")

    y = np.arange(len(sets))
    for i in x:
        ax_matrix.scatter(np.repeat(i, len(sets)), y, s=16, color="#d9dee3", zorder=1)
        active = np.flatnonzero(matrix[i])
        ax_matrix.scatter(np.repeat(i, len(active)), active, s=16, color=color, zorder=3)
        if len(active) > 1:
            ax_matrix.plot([i, i], [active.min(), active.max()], color=color, linewidth=0.9, zorder=2)
    ax_matrix.set_ylim(len(sets) - 0.6, -0.6)
    ax_matrix.set_yticks(np.arange(len(sets)))
    ax_matrix.set_yticklabels([display(x) for x in sets], fontsize=5.0)
    ax_matrix.set_xticks(x)
    ax_matrix.set_xticklabels([str(i + 1) for i in x], fontsize=4.8)
    ax_matrix.set_xlabel("Intersection", labelpad=2)
    ax_matrix.set_ylabel(set_label, labelpad=3)
    ax_matrix.tick_params(length=0, pad=2)
    for spine in ax_matrix.spines.values():
        spine.set_visible(False)


def save(fig, path: Path):
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--merged", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    tissue_gene_tissues, tissue_gene_cells = load_membership(args.merged)
    across = build_across(tissue_gene_tissues)
    within = build_within(tissue_gene_cells)
    across.to_csv(args.out_dir / "Source_01_across_tissue_ieGene_upset.tsv", sep="\t", index=False)
    within.to_csv(args.out_dir / "Source_02_within_tissue_ieGene_upset.tsv", sep="\t", index=False)

    multi = across[across["n_tissues"] >= 2].copy()
    fig = plt.figure(figsize=(8.8, 5.5))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.25], hspace=0.42)
    top = gs[0].subgridspec(2, 1, height_ratios=[1.0, 1.15], hspace=0.02)
    draw_upset(
        fig.add_subplot(top[0]), fig.add_subplot(top[1]), multi,
        "tissue_set", "Across-tissue ieGene sharing",
        "Same ieGenes observed in at least two tissues; variants collapsed within gene",
        "Tissue", color="#496b87",
    )
    fig.text(0.015, 0.985, "a", fontsize=10, fontweight="bold", va="top")
    save(fig, args.out_dir / "01_across_tissue_ieGene_upset")

    informative = []
    for tissue in sorted(within["tissue"].unique()):
        group = within[within["tissue"] == tissue]
        if group["n_shared_ieGenes"].sum() > 0:
            informative.append(tissue)
    n = len(informative)
    cols = 2
    rows_n = int(np.ceil(n / cols))
    fig = plt.figure(figsize=(9.0, max(4.8, rows_n * 3.1)))
    outer = fig.add_gridspec(rows_n, cols, hspace=0.65, wspace=0.32)
    for idx, tissue in enumerate(informative):
        r, c = divmod(idx, cols)
        cell = outer[r, c].subgridspec(2, 1, height_ratios=[0.8, 1.0], hspace=0.02)
        group = within[within["tissue"] == tissue]
        draw_upset(
            fig.add_subplot(cell[0]), fig.add_subplot(cell[1]), group,
            "cell_type_set", display(tissue),
            "Same ieGenes detected in multiple cell types",
            "Cell type", color=TISSUE_COLORS.get(display(tissue), "#496b87"),
        )
    fig.text(0.015, 0.995, "b", fontsize=10, fontweight="bold", va="top")
    fig.text(0.06, 0.02, "Only tissues with at least one within-tissue shared ieGene are shown.", fontsize=6.2, color="#52606d")
    save(fig, args.out_dir / "02_within_tissue_ieGene_upset")
    print(f"across_multi_tissue_intersections={len(multi)}")
    print(f"within_informative_tissues={len(informative)}")
    print("within_tissues=" + ",".join(display(x) for x in informative))


if __name__ == "__main__":
    main()
