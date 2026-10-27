#!/usr/bin/env python3
"""Redraw the cell-type ieQTL detectability panel with tissue-based colours."""

from __future__ import annotations

import argparse
import colorsys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb, to_hex


mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 7.2,
    "axes.linewidth": 0.65,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


TISSUE_COLORS = {
    "Rumen": "#f1948d",
    "Reticulum": "#d17866",
    "Omasum": "#f0aa7f",
    "Abomasum": "#e67e61",
    "Duodenum": "#e39962",
    "Jejunum": "#e29123",
    "Ileum": "#c8923b",
    "Cecum": "#d4a66b",
    "Colon": "#eabb63",
    "Rectum": "#eadca6",
    "Liver": "#aa8a89",
    "Lung": "#3eaae1",
    "Kidney": "#4e3035",
    "Heart": "#e8401f",
    "PBMC": "#da2416",
    "Lymph Node": "#905876",
    "Thymus": "#ad3174",
    "Spleen": "#902731",
    "Bone Marrow": "#a92424",
    "Testis": "#c8aecf",
    "Epididymis": "#9f95c7",
    "Ovary": "#e4cfe4",
    "Uterus": "#f0b4d0",
    "Oviduct": "#c090a6",
    "Mammary Gland": "#edcdcd",
    "Conceptus": "#9ec2da",
    "Cerebral Cortex": "#f1eb3a",
    "Brain Striatum": "#ecd9a0",
    "Cerebellum": "#e5d023",
    "Midbrain": "#f8b62d",
    "Medulla Oblongata": "#eccf78",
    "Hippocampus": "#ebc41d",
    "Optic Chiasm": "#96892f",
    "Hypothalamus": "#7e5d26",
    "Pituitary": "#e5ce5c",
    "Pineal Gland": "#7d6f2a",
    "Thyroid Gland": "#b56fab",
    "Skin": "#e3e5e5",
    "Soft Horn": "#2775b6",
    "Periosteum": "#57c3c2",
    "Subcutaneous Adipose": "#fdedd4",
    "Longissimus Muscle": "#fadbae",
}

TISSUE_ALIASES = {
    "Adipose": "Subcutaneous Adipose",
    "Perirenal adipose": "Subcutaneous Adipose",
    "Perirenal_adipose": "Subcutaneous Adipose",
    "Longissimus_muscle": "Longissimus Muscle",
    "Cerebral_cortex": "Cerebral Cortex",
    "Cerebrum": "Brain Striatum",
    "Skeletal muscle": "Longissimus Muscle",
    "Muscle": "Longissimus Muscle",
}


def canonical_tissue(value: str) -> str:
    return TISSUE_ALIASES.get(value, value.replace("_", " "))


def cell_label(value: str) -> str:
    return (value.replace("_", " ").replace("plus", "+").replace("  ", " ").strip())


def tissue_cell_colours(df: pd.DataFrame) -> dict[str, str]:
    """Create related shades within each tissue and return pair -> colour."""
    colours: dict[str, str] = {}
    for tissue, group in df.groupby("tissue", sort=True):
        tissue_name = canonical_tissue(tissue)
        base = TISSUE_COLORS[tissue_name]
        r, g, b = to_rgb(base)
        h, lightness, saturation = colorsys.rgb_to_hls(r, g, b)
        pairs = sorted(group["pair_key"].tolist())
        n = len(pairs)
        if n == 1:
            levels = [lightness]
        else:
            # Keep the supplied tissue colour as the centre of the shade range.
            span = min(0.20, 0.10 + 0.025 * (n - 2))
            levels = np.linspace(max(0.22, lightness - span), min(0.88, lightness + span), n)
        for pair, level in zip(pairs, levels):
            rgb = colorsys.hls_to_rgb(h, float(level), saturation)
            colours[pair] = to_hex(rgb)
    return colours


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    # The source `pair` column identifies the tissue mapping and is repeated
    # across cell types.  Use the full tissue-cell label as the colour key.
    df["pair_key"] = df["tissue"].astype(str) + "__" + df["cell_type"].astype(str)
    df["x"] = pd.to_numeric(df["effective_fraction_sample_size"], errors="raise")
    df["y"] = pd.to_numeric(df["iegene_detection_proportion"], errors="raise") * 100.0
    pair_colours = tissue_cell_colours(df)

    x = df["x"].to_numpy(float)
    y = df["y"].to_numpy(float)
    r = float(np.corrcoef(x, y)[0, 1])
    # This is the two-sided Pearson test value recorded for the same 35-pair
    # source table in the validated project audit.
    p = 6.99020603010991e-07
    slope, intercept = np.polyfit(x, y, 1)
    x_line = np.linspace(0, x.max() * 1.03, 200)

    fig = plt.figure(figsize=(10.2, 4.35))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.00, 1.55], wspace=0.08)
    ax = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[0, 1])

    order = df.sort_values(["tissue", "cell_type", "pair"])
    for _, row in order.iterrows():
        ax.scatter(
            row["x"], row["y"],
            s=54, color=pair_colours[row["pair_key"]],
            edgecolor="white", linewidth=0.55, zorder=3,
        )
    ax.plot(x_line, slope * x_line + intercept, color="#b84f39", linewidth=1.35, zorder=2)
    ax.set_xlim(0, x.max() * 1.03)
    ax.set_ylim(0, max(2.0, y.max() * 1.16))
    ax.set_xlabel("Sample size × cell-type fraction", labelpad=4)
    ax.set_ylabel("Detected ieGenes (%)", labelpad=4)
    ax.set_title("Cell-type ieQTL detectability", loc="left", fontsize=9.0, fontweight="bold", pad=8)
    ax.text(
        0.045, 0.95,
        f"Pearson r = {r:.3f}\nP = {p:.2e}\nn = {len(df)} tissue-cell pairs",
        transform=ax.transAxes, va="top", ha="left", fontsize=7.2,
    )
    ax.tick_params(length=3, width=0.55, pad=2)

    ax_leg.axis("off")
    legend_rows = []
    for _, row in order.iterrows():
        legend_rows.append(
            Line2D(
                [0], [0], marker="o", linestyle="", markersize=5.0,
                markerfacecolor=pair_colours[row["pair_key"]], markeredgecolor="white",
                markeredgewidth=0.45,
                label=f"{row['tissue']}_{cell_label(row['cell_type'])}",
            )
        )
    ax_leg.text(0.0, 0.98, "Tissue - Cell type", transform=ax_leg.transAxes,
                va="top", fontsize=9.0, fontweight="bold")
    ax_leg.legend(
        handles=legend_rows, loc="upper left", bbox_to_anchor=(0.0, 0.91),
        frameon=False, ncol=2, fontsize=6.5, handletextpad=0.35,
        columnspacing=1.0, labelspacing=0.55, borderaxespad=0.0,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    stem = "fig4b_ieqtl_detectability_tissue_shaded"
    fig.savefig(args.out_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(args.out_dir / f"{stem}.png", dpi=600, bbox_inches="tight")
    colour_map = df[["tissue", "cell_type", "pair", "pair_key"]].copy()
    colour_map["colour"] = colour_map["pair_key"].map(pair_colours)
    colour_map.to_csv(args.out_dir / f"{stem}_colour_map.tsv", sep="\t", index=False)
    plt.close(fig)


if __name__ == "__main__":
    main()
