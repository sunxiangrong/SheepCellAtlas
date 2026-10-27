#!/usr/bin/env python3
"""Reference-style lineage-stratified sheep cell type--trait dot plots.

Every mapped tissue--cell type result for the selected traits is drawn. Points
are horizontally dispersed within each lineage rather than collapsed at a
single x coordinate. The y-axis is -log10(merged_BG_adjp), with the FDR 0.05
threshold shown explicitly. Pre-specified contexts from the main/supplementary
examples are annotated; annotations do not select the plotted population.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


TRAITS = ["ALLGE3", "P3R_rint", "OHW", "CD"]
TRAIT_LABELS = {
    "ALLGE3": "ALLGE3 | reproduction",
    "P3R_rint": "P3R_rint | reproduction",
    "OHW": "OHW | wool",
    "CD": "CD | body size",
}
LINEAGES = [
    "Endothelial cells",
    "Epithelial cells",
    "Germline cells",
    "Immune cells",
    "Muscle cells",
    "Nervous cells",
    "Other cells",
    "Stromal cells",
]
LINEAGE_COLORS = {
    "Endothelial cells": "#4C9E52",
    "Epithelial cells": "#C94F7C",
    "Germline cells": "#2F6B4F",
    "Immune cells": "#C77E45",
    "Muscle cells": "#B795C9",
    "Nervous cells": "#E6C229",
    "Other cells": "#9A9A9A",
    "Stromal cells": "#6FB3D2",
}
FDR_THRESHOLD = 0.05
ANNOTATIONS = {
    "ALLGE3": [
        ("Conceptus", "Hypoblasts cells", "Hypoblasts\nConceptus", 0.24, 0.42),
        ("Oviduct", "Luminal cells", "Luminal cells\nOviduct", -0.30, 0.32),
    ],
    "P3R_rint": [
        ("Uterus", "γδT cells", "gamma-delta T cells\nUterus", -0.38, 0.27),
        ("Cerebralcortex", "SST+ GABAergic neurons", "SST+ GABA neurons\nCerebral cortex", 0.25, 0.34),
    ],
    "OHW": [
        ("PBMC", "Non-classical monocytes", "Non-classical monocytes\nPBMC", 0.22, 0.30),
    ],
    "CD": [
        ("Lung", "CD4+ T cells", "CD4+ T cells\nLung", -0.30, 0.26),
        ("Brainstriatum", "Oligodendrocyte precursor cells", "OPCs\nBrain striatum", 0.25, 0.34),
    ],
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 8,
            "axes.linewidth": 0.75,
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def read_source(result_root: Path, lineage_map_path: Path, traits: list[str]) -> pd.DataFrame:
    """Read all selected tissue--trait summaries from portable user inputs."""
    lineage = pd.read_csv(lineage_map_path, sep="\t")
    lineage_map = dict(zip(lineage["cell_type"].astype(str), lineage["cell_lineage"].astype(str)))
    tables: list[pd.DataFrame] = []
    for path in result_root.glob("*/*/Merged_celltype_summary.csv"):
        tissue, trait = path.parent.parent.name, path.parent.name
        if trait not in traits:
            continue
        table = pd.read_csv(path, usecols=["celltype", "mean_TRS", "merged_BG_adjp", "n_cells"])
        table["tissue"] = tissue
        table["trait"] = trait
        table["source_summary"] = str(path)
        tables.append(table)
    out = pd.concat(tables, ignore_index=True)
    out["mean_TRS"] = pd.to_numeric(out["mean_TRS"], errors="coerce")
    out["merged_BG_adjp"] = pd.to_numeric(out["merged_BG_adjp"], errors="coerce")
    out["lineage"] = out["celltype"].map(lineage_map).fillna("Unmapped")
    out = out.loc[out["lineage"].isin(LINEAGES) & out["merged_BG_adjp"].notna()].copy()
    out["neglog10_fdr"] = -np.log10(out["merged_BG_adjp"].clip(lower=np.nextafter(0, 1)))
    out["fdr_significant"] = out["merged_BG_adjp"] < FDR_THRESHOLD
    out["celltype_context"] = out["tissue"] + " | " + out["celltype"]
    return out


def add_dispersed_x(frame: pd.DataFrame, seed: int = 20260921) -> pd.DataFrame:
    """Give all points a reproducible, compact horizontal spread per lineage."""
    rng = np.random.default_rng(seed)
    pieces = []
    lineage_index = {name: i for i, name in enumerate(LINEAGES)}
    for (trait, lineage), group in frame.groupby(["trait", "lineage"], sort=False):
        group = group.sort_values(["tissue", "celltype"], kind="stable").copy()
        n = len(group)
        # Offsets are intentionally independent of FDR. This makes a compact
        # strip/beeswarm-like cloud without creating a false diagonal trend.
        offsets = rng.uniform(-0.31, 0.31, n) if n > 1 else np.array([0.0])
        group["x"] = lineage_index[lineage] + offsets
        pieces.append(group)
    return pd.concat(pieces, ignore_index=True)


def annotate_examples(ax: plt.Axes, frame: pd.DataFrame, trait: str) -> None:
    for tissue, celltype, label, dx, dy in ANNOTATIONS[trait]:
        row = frame.loc[(frame["tissue"].eq(tissue)) & (frame["celltype"].eq(celltype))]
        if row.empty:
            continue
        point = row.sort_values("neglog10_fdr", ascending=False).iloc[0]
        color = LINEAGE_COLORS[point["lineage"]]
        ax.annotate(
            label,
            xy=(point["x"], point["neglog10_fdr"]),
            xytext=(point["x"] + dx, point["neglog10_fdr"] + dy),
            fontsize=6.5,
            color=color,
            ha="left" if dx >= 0 else "right",
            va="bottom",
            arrowprops={"arrowstyle": "-", "lw": 0.55, "color": color, "shrinkA": 0, "shrinkB": 2},
            zorder=5,
        )


def draw(data: pd.DataFrame, output: Path, traits: list[str]) -> None:
    data = add_dispersed_x(data)
    fig, axes = plt.subplots(len(traits), 1, figsize=(7.45, 3.25 * len(traits)), constrained_layout=False)
    if len(traits) == 1:
        axes = [axes]
    fig.subplots_adjust(left=0.12, right=0.80, top=0.965, bottom=0.065, hspace=0.42)
    threshold_y = -np.log10(FDR_THRESHOLD)
    for panel_i, (ax, trait) in enumerate(zip(axes, traits)):
        frame = data.loc[data["trait"].eq(trait)].copy()
        # All contexts retain their lineage identity. Non-significant points are
        # lighter and smaller, rather than being omitted or recoloured gray.
        nonsig = frame.loc[~frame["fdr_significant"]]
        sig = frame.loc[frame["fdr_significant"]]
        for lineage in LINEAGES:
            rows = nonsig.loc[nonsig["lineage"].eq(lineage)]
            ax.scatter(rows["x"], rows["neglog10_fdr"], s=12, c=LINEAGE_COLORS[lineage], alpha=0.30, linewidths=0, zorder=1, rasterized=True)
        for lineage in LINEAGES:
            rows = sig.loc[sig["lineage"].eq(lineage)]
            ax.scatter(rows["x"], rows["neglog10_fdr"], s=19, c=LINEAGE_COLORS[lineage], alpha=0.90, edgecolors="none", zorder=2, rasterized=True)
        ax.axhline(threshold_y, color="#353535", lw=0.75, ls=(0, (4, 2.5)), zorder=0)
        ax.text(len(LINEAGES) - 0.42, threshold_y + 0.10, "FDR = 0.05", ha="right", va="bottom", fontsize=6.3, color="#333333")
        annotate_examples(ax, frame, trait)
        ymax = max(float(frame["neglog10_fdr"].max()) + 0.75, threshold_y + 1.15)
        ax.set_ylim(-0.12, ymax)
        ax.set_xlim(-0.58, len(LINEAGES) - 0.42)
        ax.set_xticks(range(len(LINEAGES)))
        ax.set_xticklabels([label.replace(" cells", "") for label in LINEAGES], rotation=40, ha="right", fontsize=6.7)
        ax.tick_params(axis="x", length=0, pad=2)
        ax.tick_params(axis="y", labelsize=7, width=0.65, length=3)
        ax.set_ylabel(r"$-\log_{10}$(FDR)", fontsize=8)
        ax.set_title(TRAIT_LABELS[trait], fontsize=9.5, fontweight="bold", pad=5)
        ax.text(-0.10, 1.03, chr(ord("a") + panel_i), transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_linewidth(0.75)
        ax.spines["left"].set_linewidth(0.75)
        if panel_i == len(traits) - 1:
            ax.set_xlabel("Cell lineage", fontsize=8.2)
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=LINEAGE_COLORS[name], markersize=4.6, label=name)
        for name in LINEAGES
    ]
    handles.extend([
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#555555", alpha=0.9, markersize=5.0, label="FDR < 0.05"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#555555", alpha=0.30, markersize=4.2, label="FDR ≥ 0.05"),
    ])
    fig.legend(handles=handles, title="Cell lineage", title_fontsize=7.5, fontsize=7.0, frameon=False, loc="center left", bbox_to_anchor=(0.815, 0.52), borderaxespad=0.0, labelspacing=0.75)
    fig.savefig(output.with_suffix(".png"), dpi=500, bbox_inches="tight", facecolor="white")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-root", type=Path, required=True,
                        help="Directory containing <tissue>/<trait>/Merged_celltype_summary.csv files")
    parser.add_argument("--lineage-map", type=Path, required=True,
                        help="TSV with cell_type and cell_lineage columns")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--traits", nargs="+", default=TRAITS, choices=TRAITS)
    args = parser.parse_args()
    setup_style()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    data = read_source(args.result_root, args.lineage_map, args.traits)
    draw(data, args.out_dir / "sheep_lineage_celltype_fdr_reference_style", args.traits)
    data.sort_values(["trait", "lineage", "tissue", "celltype"]).to_csv(args.out_dir / "source_all_contexts.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "traits": ";".join(args.traits),
        "point_unit": "one tissue-cell type-trait association from Merged_celltype_summary.csv",
        "x_axis": "mapped cell lineage in fixed alphabetical order",
        "y_axis": "-log10(merged_BG_adjp)",
        "significance_line": f"FDR={FDR_THRESHOLD}; y={-np.log10(FDR_THRESHOLD):.12g}",
        "point_colours": "lineage for all contexts; FDR >= 0.05 points use lower opacity and smaller size",
        "horizontal_spread": "deterministic within-trait-lineage offsets, seed 20260921",
        "labels": "pre-specified main/supplementary example contexts only",
        "n_plotted_contexts": len(data),
        "n_fdr_significant_contexts": int(data["fdr_significant"].sum()),
    }]).to_csv(args.out_dir / "run_metadata.tsv", sep="\t", index=False)
    print(f"plotted={len(data)} significant={int(data['fdr_significant'].sum())} out={args.out_dir}")


if __name__ == "__main__":
    main()
