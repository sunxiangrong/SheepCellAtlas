#!/usr/bin/env python3
"""Create preliminary, matched GWAS--cell-type ieQTL colocalization locus plots."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
        "axes.spines.right": False,
        "axes.spines.top": False,
    }
)

GWAS_COLOR = "#3F6F9F"
IEQTL_COLOR = "#C75B5B"
BACKGROUND = "#B8BEC5"
LEAD_COLOR = "#222222"


def neglog10(values: pd.Series) -> pd.Series:
    return -np.log10(values.clip(lower=np.finfo(float).tiny))


def plot_candidate(row: pd.Series, args: argparse.Namespace) -> None:
    chrom, lead_pos = row.lead_snp.split("_", 1)
    lead_pos = int(lead_pos)
    start = max(0, lead_pos - args.window_bp)
    end = lead_pos + args.window_bp

    gwas_file = Path(args.gwas_root) / row.trait / f"filtered_STDERR_{row.trait}1.tbl.txt"
    ieqtl_file = Path(args.ieqtl_root) / row.tissue / f"{row.cell_type}.cis_qtl.tsv.gz"
    if not gwas_file.is_file() or not ieqtl_file.is_file():
        raise FileNotFoundError(f"Missing input for {row.trait}, {row.tissue}, {row.cell_type}")

    gwas = pd.read_csv(gwas_file, sep="\t", usecols=["Chr", "SNP", "bp", "p"])
    gwas = gwas[(gwas["Chr"].astype(str) == str(chrom)) & gwas["bp"].between(start, end)].copy()
    gwas["logp"] = neglog10(gwas["p"])

    ieqtl = pd.read_csv(
        ieqtl_file,
        sep="\t",
        compression="gzip",
        usecols=["pheno_id", "variant_id", "beta_g2", "beta_se_g2", "pval_g2"],
    )
    ieqtl = ieqtl[ieqtl["pheno_id"] == row.gene].copy()
    coords = ieqtl["variant_id"].str.split("_", n=1, expand=True)
    ieqtl["chrom"] = coords[0]
    ieqtl["bp"] = pd.to_numeric(coords[1], errors="coerce")
    ieqtl = ieqtl[(ieqtl["chrom"] == str(chrom)) & ieqtl["bp"].between(start, end)].copy()
    ieqtl["logp"] = neglog10(ieqtl["pval_g2"])

    if gwas.empty or ieqtl.empty:
        raise ValueError(f"No locus data retained for {row.trait}, {row.tissue}, {row.cell_type}, {row.gene}")

    fig, axes = plt.subplots(2, 1, figsize=(4.25, 4.1), sharex=True, gridspec_kw={"hspace": 0.13})
    gwas_lead = gwas[gwas["SNP"] == row.lead_snp]
    if gwas_lead.empty:
        gwas_lead = gwas.loc[[gwas["p"].idxmin()]]
    ieqtl_lead = ieqtl.loc[[ieqtl["pval_g2"].idxmin()]]
    gwas_lead_pos = float(gwas_lead["bp"].iloc[0])
    ieqtl_lead_pos = float(ieqtl_lead["bp"].iloc[0])

    plot_specs = [
        (axes[0], gwas, "bp", "logp", "GWAS", GWAS_COLOR, gwas_lead),
        (axes[1], ieqtl, "bp", "logp", "Cell-type ieQTL interaction", IEQTL_COLOR, ieqtl_lead),
    ]
    for ax, data, x, y, label, color, lead in plot_specs:
        ax.scatter(data[x] / 1e6, data[y], color=BACKGROUND, s=10, linewidths=0, alpha=0.82, rasterized=True)
        ax.scatter(lead[x] / 1e6, lead[y], color=color, edgecolor=LEAD_COLOR, linewidth=0.45,
                   marker="D", s=32, zorder=3)
        ax.axvline(gwas_lead_pos / 1e6, color=GWAS_COLOR, linewidth=0.8, linestyle="--", alpha=0.75)
        ax.axvline(ieqtl_lead_pos / 1e6, color=IEQTL_COLOR, linewidth=0.8, linestyle=":", alpha=0.85)
        ax.set_ylabel(r"$-\log_{10}(P)$")
        ax.set_title(label, loc="left", fontsize=8, fontweight="bold", pad=3)
        ax.tick_params(direction="out", length=3, width=0.7)

    axes[1].set_xlabel(f"Position on chr{chrom} (Mb)")
    axes[0].text(
        0.99,
        0.95,
        f"{row.trait} | {row.tissue} | {row.cell_type}\n{row.gene}\n"
        + r"PP$_{H4}$" + f" = {row.pph4:.3f}; " + r"PP$_{H3}$" + f" = {row.pph3:.3f}",
        transform=axes[0].transAxes,
        ha="right",
        va="top",
        fontsize=7.2,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.92, "pad": 2.3},
    )
    axes[1].text(
        0.01,
        0.93,
        r"interaction effect ($\beta_{g2}$)",
        transform=axes[1].transAxes,
        ha="left",
        va="top",
        fontsize=7.2,
        color="#555555",
    )
    axes[1].text(
        0.99,
        0.93,
        f"GWAS lead: {gwas_lead['SNP'].iloc[0]}\n"
        f"ieQTL lead: {ieqtl_lead['variant_id'].iloc[0]}",
        transform=axes[1].transAxes,
        ha="right",
        va="top",
        fontsize=6.8,
        color="#444444",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 1.8},
    )
    fig.align_ylabels(axes)

    stem = f"{row.trait}__{row.tissue}__{row.cell_type}__{row.gene}__{row.lead_snp}"
    stem = stem.replace("/", "_").replace(" ", "_")
    outbase = Path(args.out_dir) / stem
    fig.savefig(outbase.with_suffix(".png"), dpi=350, bbox_inches="tight", facecolor="white")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--gwas-root", required=True)
    parser.add_argument("--ieqtl-root", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--window-bp", type=int, default=500_000)
    args = parser.parse_args()

    candidates = pd.read_csv(args.candidates, sep="\t")
    required = {"trait", "tissue", "cell_type", "gene", "lead_snp", "pph3", "pph4"}
    missing = required.difference(candidates.columns)
    if missing:
        raise ValueError(f"Candidate file lacks columns: {sorted(missing)}")
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    for _, row in candidates.iterrows():
        plot_candidate(row, args)


if __name__ == "__main__":
    main()
