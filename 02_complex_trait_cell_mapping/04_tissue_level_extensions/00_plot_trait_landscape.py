#!/usr/bin/env python3
"""Draw a reference-style sheep trait tissue-cell type dot plot with variable lineage widths.

Each point is one tissue--cell type association from the current tissue-level
scPagwas summaries. Lineage blocks are allotted horizontal width according to
their observed point count, so dense groups are not forced into the same space
as sparse groups such as Germline cells.
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


FDR_THRESHOLD = 0.05
LINEAGES = [
    "Endothelial cells", "Epithelial cells", "Germline cells", "Immune cells",
    "Muscle cells", "Nervous cells", "Other cells", "Stromal cells",
]
LINEAGE_COLORS = {
    "Endothelial cells": "#14824B",
    "Epithelial cells": "#D94C98",
    "Germline cells": "#39A935",
    "Immune cells": "#F2A64A",
    "Muscle cells": "#D9A8CF",
    "Nervous cells": "#F2D319",
    "Other cells": "#9B9B9B",
    "Stromal cells": "#58BED8",
}
TRAIT_LABELS = {
    "ALLGE3": "ALLGE3 | reproduction",
    "P3R_rint": "P3R_rint | reproduction",
    "OHW": "OHW | wool",
    "CD": "CD | body size",
}
ANNOTATIONS = {
    "ALLGE3": [
        ("Conceptus", "Hypoblasts cells", "Hypoblasts\n(Conceptus)", 0.37, 0.52),
        ("Conceptus", "Trophoblast progenitors", "Trophoblast progenitors\n(Conceptus)", 0.56, 0.78),
        ("Conceptus", "Trophectoderm cells", "Trophectoderm\n(Conceptus)", -0.58, 0.10),
        ("Oviduct", "Luminal cells", "Luminal cells\n(Oviduct)", 0.34, -0.63),
    ],
    "P3R_rint": [
        ("Uterus", "γδT cells", "gamma-delta T cells\n(Uterus)", -0.43, 0.38),
        ("Uterus", "CD8+ T cells", "CD8+ T cells\n(Uterus)", 0.38, 0.20),
        ("Uterus", "NK T cells", "NK T cells\n(Uterus)", -0.34, -0.56),
        ("Uterus", "B cells", "B cells\n(Uterus)", 0.38, -0.58),
    ],
    "OHW": [
        ("PBMC", "Non-classical monocytes", "Non-classical monocytes\n(PBMC)", 0.48, 0.44),
        ("PBMC", "CD4+ T cells", "CD4+ T cells\n(PBMC)", -0.42, 0.26),
    ],
    "CD": [
        ("Lung", "CD4+ T cells", "CD4+ T cells\n(Lung)", -0.38, 0.35),
        ("Brainstriatum", "Oligodendrocyte precursor cells", "OPCs\n(Brain striatum)", 0.45, 0.36),
    ],
}


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 8,
        "axes.linewidth": 0.8,
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })


def read_data(result_root: Path, lineage_map_path: Path, trait: str) -> pd.DataFrame:
    """Read per-tissue scPagwas summaries from a user-supplied result root."""
    mapping = pd.read_csv(lineage_map_path, sep="\t")
    lineage_map = dict(zip(mapping["cell_type"].astype(str), mapping["cell_lineage"].astype(str)))
    tables = []
    for path in result_root.glob(f"*/{trait}/Merged_celltype_summary.csv"):
        tissue = path.parent.parent.name
        table = pd.read_csv(path, usecols=["celltype", "mean_TRS", "merged_BG_adjp", "n_cells"])
        table["tissue"] = tissue
        table["source_summary"] = str(path)
        tables.append(table)
    out = pd.concat(tables, ignore_index=True)
    out["mean_TRS"] = pd.to_numeric(out["mean_TRS"], errors="coerce")
    out["merged_BG_adjp"] = pd.to_numeric(out["merged_BG_adjp"], errors="coerce")
    out["lineage"] = out["celltype"].map(lineage_map).fillna("Unmapped")
    out = out.loc[out["lineage"].isin(LINEAGES) & out["merged_BG_adjp"].notna()].copy()
    out["neglog10_fdr"] = -np.log10(out["merged_BG_adjp"].clip(lower=np.nextafter(0, 1)))
    out["fdr_significant"] = out["merged_BG_adjp"] < FDR_THRESHOLD
    return out


def assign_variable_width_x(data: pd.DataFrame, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = data.groupby("lineage").size().reindex(LINEAGES, fill_value=0).rename("n_points").reset_index()
    # Square-root scaling gives dense lineages visibly more room without letting
    # one large lineage dominate the full figure width.
    counts["block_width"] = 0.85 + 0.18 * np.sqrt(counts["n_points"])
    counts["block_start"] = counts["block_width"].shift(fill_value=0).cumsum() + np.arange(len(counts)) * 0.22
    counts["block_center"] = counts["block_start"] + counts["block_width"] / 2
    lookup = counts.set_index("lineage")
    rng = np.random.default_rng(seed)
    pieces = []
    for lineage in LINEAGES:
        rows = data.loc[data["lineage"].eq(lineage)].copy()
        if rows.empty:
            continue
        start = lookup.loc[lineage, "block_start"]
        width = lookup.loc[lineage, "block_width"]
        rows["x"] = start + rng.uniform(0.06 * width, 0.94 * width, len(rows))
        pieces.append(rows)
    return pd.concat(pieces, ignore_index=True), counts


def annotate_examples(ax: plt.Axes, data: pd.DataFrame, trait: str) -> None:
    for tissue, celltype, label, dx, dy in ANNOTATIONS.get(trait, []):
        row = data.loc[data["tissue"].eq(tissue) & data["celltype"].eq(celltype)]
        if row.empty:
            continue
        point = row.iloc[0]
        color = LINEAGE_COLORS[point["lineage"]]
        ax.annotate(
            label,
            xy=(point["x"], point["neglog10_fdr"]),
            xytext=(point["x"] + dx, point["neglog10_fdr"] + dy),
            ha="left" if dx > 0 else "right",
            va="bottom",
            color=color,
            fontsize=7.2,
            arrowprops={"arrowstyle": "->", "lw": 0.65, "color": color, "shrinkA": 2, "shrinkB": 2},
            zorder=6,
        )


def scaled_point_sizes(values: pd.Series, reference: float) -> np.ndarray:
    """Encode positive raw mean_TRS without treating TRS as a probability."""
    positive = pd.to_numeric(values, errors="coerce").fillna(0).clip(lower=0)
    return 13 + 46 * np.clip(positive / reference, 0, 1)


def draw(data: pd.DataFrame, blocks: pd.DataFrame, output: Path, trait: str) -> None:
    fig, ax = plt.subplots(figsize=(10.4, 4.9))
    fig.subplots_adjust(left=0.10, right=0.70, bottom=0.16, top=0.88)
    nonsig = data.loc[~data["fdr_significant"]]
    sig = data.loc[data["fdr_significant"]]
    trs_reference = max(float(data["mean_TRS"].clip(lower=0).quantile(0.95)), 0.6)
    for lineage in LINEAGES:
        color = LINEAGE_COLORS[lineage]
        low = nonsig.loc[nonsig["lineage"].eq(lineage)]
        high = sig.loc[sig["lineage"].eq(lineage)]
        ax.scatter(low["x"], low["neglog10_fdr"], s=scaled_point_sizes(low["mean_TRS"], trs_reference), color=color, alpha=0.28, edgecolors="none", zorder=1, rasterized=True)
        ax.scatter(high["x"], high["neglog10_fdr"], s=scaled_point_sizes(high["mean_TRS"], trs_reference), color=color, alpha=0.93, edgecolors="none", zorder=2, rasterized=True)
    threshold_y = -np.log10(FDR_THRESHOLD)
    ax.axhline(threshold_y, color="#181818", lw=1.6, ls=(0, (5, 4)), zorder=0)
    ax.text(blocks["block_start"].iloc[-1] + blocks["block_width"].iloc[-1], threshold_y + 0.08, "FDR = 0.05", ha="right", va="bottom", fontsize=7.0)
    annotate_examples(ax, data, trait)
    xmax = (blocks["block_start"] + blocks["block_width"]).max()
    ax.set_xlim(-0.15, xmax + 0.15)
    ax.set_ylim(-0.12, max(float(data["neglog10_fdr"].max()) + 0.65, threshold_y + 1.0))
    ax.set_xticks(blocks["block_center"])
    ax.set_xticklabels([])
    ax.tick_params(axis="x", length=4, width=0.75)
    ax.tick_params(axis="y", labelsize=8, length=3.5, width=0.75)
    ax.set_ylabel(r"$-\log_{10}$(FDR)", fontsize=9)
    ax.set_xlabel("Cell lineage", fontsize=9, labelpad=7)
    ax.set_title(TRAIT_LABELS.get(trait, trait), fontsize=10.5, fontweight="bold", pad=12)
    ax.text(-0.095, 1.03, "a", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=LINEAGE_COLORS[x], markersize=5.0, label=x)
        for x in LINEAGES
    ]
    lineage_legend = ax.legend(handles=handles, title="Cell lineage", title_fontsize=7.5, fontsize=7.1, frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.56), labelspacing=0.68, borderaxespad=0)
    ax.add_artist(lineage_legend)
    size_values = [0.0, 0.6, 1.0]
    size_handles = [
        ax.scatter([], [], s=scaled_point_sizes(pd.Series([value]), trs_reference)[0], color="#595959", alpha=0.9, edgecolors="none")
        for value in size_values
    ]
    ax.legend(size_handles, ["0", "0.6", "1.0"], title="Mean TRS", title_fontsize=7.5, fontsize=7.1, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 0.25), labelspacing=0.9, borderaxespad=0)
    fig.savefig(output.with_suffix(".png"), dpi=500, bbox_inches="tight", facecolor="white")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait", required=True, choices=sorted(TRAIT_LABELS))
    parser.add_argument("--result-root", type=Path, required=True,
                        help="Directory containing <tissue>/<trait>/Merged_celltype_summary.csv files")
    parser.add_argument("--lineage-map", type=Path, required=True,
                        help="TSV with cell_type and cell_lineage columns")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260925)
    args = parser.parse_args()
    setup_style()
    args.outdir.mkdir(parents=True, exist_ok=True)
    data = read_data(args.result_root, args.lineage_map, args.trait)
    data, blocks = assign_variable_width_x(data, args.seed)
    output = args.outdir / f"{args.trait}_landscape"
    draw(data, blocks, output, args.trait)
    data.sort_values(["lineage", "tissue", "celltype"]).to_csv(args.outdir / f"{args.trait}_landscape_source.tsv", sep="\t", index=False)
    blocks.to_csv(args.outdir / f"{args.trait}_landscape_lineage_layout.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "trait": args.trait,
        "point_unit": "one tissue-cell type association",
        "n_plotted_contexts": len(data),
        "n_fdr_significant": int(data["fdr_significant"].sum()),
        "x_axis": "cell lineages in fixed order; block width = 0.85 + 0.18*sqrt(n_points)",
        "y_axis": "-log10(merged_BG_adjp)",
        "threshold": f"FDR={FDR_THRESHOLD}; y={-np.log10(FDR_THRESHOLD):.12g}",
        "point_style": "lineage-colored for all points; point area scaled from raw positive mean_TRS; lower opacity for FDR >= 0.05",
        "annotations": "; ".join(f"{x[0]} | {x[1]}" for x in ANNOTATIONS.get(args.trait, [])),
        "random_seed_for_x_spread": args.seed,
    }]).to_csv(args.outdir / f"{args.trait}_landscape_metadata.tsv", sep="\t", index=False)
    print(f"trait={args.trait} plotted={len(data)} significant={int(data['fdr_significant'].sum())} out={args.outdir}")


if __name__ == "__main__":
    main()
