#!/usr/bin/env python3
"""Draw the final two-case PCC/stable-ortholog panel for the psoriasis story.

The panel uses the original ranked PCC table and Fisher enrichment summary.  It
shows the rolling fraction of cross-species stable orthologs along the PCC
ranking and a compact observed-versus-expected comparison for the top 5%.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7.4,
    "axes.linewidth": 0.65,
})


CASES = [
    ("psoriasis_lymphnode", "Psoriasis | Lymph node"),
    ("psoriasis_skin_softhorn_subcelltype", "Psoriasis | Skin/Softhorn subtypes"),
]
WINDOW = 450
TOP_FRACTION = 0.05
RED = "#C94A40"
GREY = "#B9C0C5"
DARK = "#2A2A2A"


def rolling_fraction(frame: pd.DataFrame) -> pd.DataFrame:
    flags = frame["is_stable_ortholog"].astype(float)
    fraction = flags.rolling(
        window=WINDOW,
        center=True,
        min_periods=max(25, WINDOW // 5),
    ).mean()
    n = len(frame)
    return pd.DataFrame({
        "rank": np.arange(1, n + 1),
        "rank_percent": np.arange(1, n + 1) / n * 100,
        "stable_fraction": fraction,
    })


def draw_row(ax_curve, ax_bar, ranked: pd.DataFrame, stat: pd.Series, title: str) -> pd.DataFrame:
    curve = rolling_fraction(ranked)
    baseline = ranked["is_stable_ortholog"].mean()
    top_n = int(stat["top_n"])
    top_pct = top_n / len(ranked) * 100

    ax_curve.axvspan(0, top_pct, color=RED, alpha=0.075, lw=0, zorder=0)
    ax_curve.axhline(baseline, color="#8C8C8C", lw=0.75, ls=(0, (2, 2)), zorder=1)
    ax_curve.plot(curve["rank_percent"], curve["stable_fraction"], color=RED, lw=1.55, zorder=2)
    ax_curve.axvline(top_pct, color=DARK, lw=0.7, ls=(0, (2, 2)), zorder=2)

    ymax = max(float(curve["stable_fraction"].max(skipna=True)), baseline) * 1.16
    ax_curve.set_ylim(0, max(0.22, min(0.40, ymax)))
    ax_curve.set_xlim(0, 100)
    ax_curve.set_yticks([0, 0.1, 0.2, 0.3])
    ax_curve.set_title(title, loc="left", fontsize=8.6, fontweight="bold", pad=2)
    ax_curve.text(
        top_pct + 1.2,
        ax_curve.get_ylim()[1] * 0.91,
        "top 5%",
        fontsize=6.7,
        color=DARK,
        va="top",
    )
    ax_curve.text(
        98.0,
        baseline + 0.006,
        "background",
        fontsize=6.3,
        color="#666666",
        ha="right",
        va="bottom",
    )
    ax_curve.spines[["top", "right"]].set_visible(False)
    ax_curve.tick_params(length=2.2, width=0.6, pad=2)

    expected = float(stat["expected_stable_in_top"])
    observed = float(stat["stable_in_top"])
    values = [expected, observed]
    colors = [GREY, RED]
    labels = ["Expected", "Observed"]
    bars = ax_bar.bar([0, 1], values, width=0.58, color=colors, edgecolor="#4A4A4A", linewidth=0.5)
    for bar, value in zip(bars, values):
        ax_bar.text(bar.get_x() + bar.get_width() / 2, value + max(values) * 0.035, f"{value:.0f}",
                    ha="center", va="bottom", fontsize=7.0)
    ax_bar.set_xticks([0, 1], labels, fontsize=6.4)
    ax_bar.tick_params(axis="x", length=0, pad=2)
    ax_bar.set_ylim(0, max(values) * 1.30)
    ax_bar.set_yticks([])
    ax_bar.spines[["top", "right", "left"]].set_visible(False)
    ax_bar.spines["bottom"].set_linewidth(0.6)
    ax_bar.text(0.02, 0.98, f"OR {stat['odds_ratio']:.2f}\nFDR {stat['fisher_fdr']:.1e}",
                transform=ax_bar.transAxes, ha="left", va="top", fontsize=6.9, color=DARK)
    return curve


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    ranked = pd.read_csv(args.input_dir / "candidate_pcc_with_stable_ortholog_flag.tsv.gz", sep="\t")
    summary = pd.read_csv(args.input_dir / "candidate_stable_fraction_enrichment_summary.tsv", sep="\t")

    fig = plt.figure(figsize=(7.05, 4.15), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[4.5, 1.15], hspace=0.24, wspace=0.10)
    curve_frames = []
    source_rows = []
    for row_index, (case_id, title) in enumerate(CASES):
        ax_curve = fig.add_subplot(grid[row_index, 0])
        ax_bar = fig.add_subplot(grid[row_index, 1])
        case_ranked = ranked.loc[ranked["case_id"].eq(case_id)].copy().reset_index(drop=True)
        case_stat = summary.loc[summary["case_id"].eq(case_id)].iloc[0]
        curve = draw_row(ax_curve, ax_bar, case_ranked, case_stat, title)
        curve.insert(0, "case_id", case_id)
        curve_frames.append(curve)
        source_rows.append(case_ranked)
        ax_curve.set_ylabel("Stable ortholog\nfraction")
        if row_index == 0:
            ax_curve.tick_params(labelbottom=False)
        else:
            ax_curve.set_xlabel("Genes ranked by PCC percentile")

    for ext in ("pdf", "svg", "png"):
        fig.savefig(
            args.out_dir / f"panelC_psoriasis_pcc_stable_ortholog.{ext}",
            dpi=450 if ext == "png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)

    selected_summary = summary.loc[summary["case_id"].isin([x[0] for x in CASES])].copy()
    selected_summary.to_csv(args.out_dir / "panelC_stable_ortholog_enrichment_stats.tsv", sep="\t", index=False)
    pd.concat(curve_frames, ignore_index=True).to_csv(args.out_dir / "panelC_rolling_stable_fraction.tsv", sep="\t", index=False)
    pd.concat(source_rows, ignore_index=True).to_csv(
        args.out_dir / "panelC_pcc_ranked_gene_source.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    print(args.out_dir)


if __name__ == "__main__":
    main()
