#!/usr/bin/env python3
"""Publication-ready locus plot for ear-length GWAS and skeletal-muscle-cell ieQTL."""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pandas as pd


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.75,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

LD_COLORS = ["#3F78A8", "#9DC3DE", "#F7D98B", "#E98D63", "#C84943"]
LD_BOUNDS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.001]
LD_CMAP = ListedColormap(LD_COLORS)
LD_NORM = BoundaryNorm(LD_BOUNDS, LD_CMAP.N)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gwas", required=True)
    parser.add_argument("--ieqtl", required=True)
    parser.add_argument("--bfile", required=True)
    parser.add_argument("--plink", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--lead", default="3_7603133")
    parser.add_argument("--trait-label", default="Ear length")
    parser.add_argument("--tissue-label", default="Muscle")
    parser.add_argument("--cell-type", default="Skeletal muscle cells")
    parser.add_argument("--gene", default="STXBP1")
    parser.add_argument("--pph3", type=float, default=0.00957286266610383)
    parser.add_argument("--pph4", type=float, default=0.947003006337563)
    parser.add_argument("--shared-snps", type=int, default=173)
    parser.add_argument("--window-bp", type=int, default=500_000)
    return parser.parse_args()


def ld_table(args: argparse.Namespace, chrom: str, start: int, end: int) -> pd.DataFrame:
    with tempfile.TemporaryDirectory() as tmp:
        prefix = Path(tmp) / "lead_ld"
        cmd = [
            args.plink,
            "--bfile", args.bfile,
            "--chr-set", "26",
            "--chr", chrom,
            "--from-bp", str(start),
            "--to-bp", str(end),
            "--ld-snp", args.lead,
            "--r2", "gz",
            "--ld-window", "999999",
            "--ld-window-kb", "1000",
            "--ld-window-r2", "0",
            "--out", str(prefix),
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        path = Path(f"{prefix}.ld.gz")
        ld = pd.read_csv(path, sep=r"\s+", compression="gzip", usecols=["SNP_B", "R2"])
    return ld.rename(columns={"SNP_B": "variant_id", "R2": "r2"})


def add_ld(frame: pd.DataFrame, ld: pd.DataFrame, variant_column: str, lead: str) -> pd.DataFrame:
    frame = frame.merge(ld, left_on=variant_column, right_on="variant_id", how="left")
    frame.loc[frame[variant_column] == lead, "r2"] = 1.0
    frame["r2"] = frame["r2"].fillna(0.0)
    return frame.sort_values("r2")


def draw_track(ax: plt.Axes, data: pd.DataFrame, variant_column: str, p_column: str, title: str, lead: str) -> None:
    ax.scatter(
        data["bp"] / 1e6,
        -np.log10(data[p_column].clip(lower=np.finfo(float).tiny)),
        c=data["r2"],
        cmap=LD_CMAP,
        norm=LD_NORM,
        s=13,
        linewidths=0,
        rasterized=True,
    )
    lead_row = data[data[variant_column] == lead]
    if not lead_row.empty:
        lead_y = -np.log10(float(lead_row[p_column].iloc[0]))
        lead_x = float(lead_row["bp"].iloc[0]) / 1e6
        ax.scatter(lead_x, lead_y, s=42, marker="D", c="#111111", edgecolors="white", linewidths=0.55, zorder=4)
    ax.set_title(title, loc="left", fontsize=7.4, fontweight="bold", pad=2)
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.tick_params(direction="out", length=3, width=0.7)


def main() -> None:
    args = parse_args()
    chrom, pos = args.lead.split("_", 1)
    center = int(pos)
    start, end = center - args.window_bp, center + args.window_bp

    gwas = pd.read_csv(args.gwas, sep="\t", usecols=["Chr", "SNP", "bp", "p"])
    gwas = gwas[(gwas["Chr"].astype(str) == chrom) & gwas["bp"].between(start, end)].copy()
    gwas = gwas.rename(columns={"SNP": "variant_id"})

    ieqtl = pd.read_csv(
        args.ieqtl,
        sep="\t",
        compression="gzip",
        usecols=["pheno_id", "variant_id", "pval_g2", "beta_g2"],
    )
    ieqtl = ieqtl[ieqtl["pheno_id"] == args.gene].copy()
    coords = ieqtl["variant_id"].str.split("_", n=1, expand=True)
    ieqtl["chrom"] = coords[0]
    ieqtl["bp"] = pd.to_numeric(coords[1], errors="coerce")
    ieqtl = ieqtl[(ieqtl["chrom"] == chrom) & ieqtl["bp"].between(start, end)].copy()
    if gwas.empty or ieqtl.empty:
        raise RuntimeError("No data in requested locus window")

    ld = ld_table(args, chrom, start, end)
    gwas = add_ld(gwas, ld, "variant_id", args.lead)
    ieqtl = add_ld(ieqtl, ld, "variant_id", args.lead)
    ieqtl_lead = ieqtl.loc[ieqtl["pval_g2"].idxmin()]

    fig = plt.figure(figsize=(6.55, 4.5), constrained_layout=False)
    grid = fig.add_gridspec(
        nrows=2,
        ncols=2,
        width_ratios=[1, 0.045],
        height_ratios=[1, 1],
        left=0.12,
        right=0.92,
        top=0.83,
        bottom=0.14,
        hspace=0.16,
        wspace=0.24,
    )
    ax_gwas = fig.add_subplot(grid[0, 0])
    ax_ieqtl = fig.add_subplot(grid[1, 0], sharex=ax_gwas)
    cax = fig.add_subplot(grid[:, 1])

    draw_track(ax_gwas, gwas, "variant_id", "p", f"{args.trait_label} GWAS", args.lead)
    draw_track(ax_ieqtl, ieqtl, "variant_id", "pval_g2", f"{args.cell_type} ieQTL interaction", args.lead)
    ax_ieqtl.set_xlabel(f"Position on chr{chrom} (Mb)")
    plt.setp(ax_gwas.get_xticklabels(), visible=False)

    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=LD_NORM, cmap=LD_CMAP), cax=cax, ticks=LD_BOUNDS[:-1])
    cb.set_label(r"LD ($r^2$)", labelpad=3)
    cb.ax.tick_params(labelsize=6.5, length=2)

    fig.text(0.12, 0.955, f"{args.trait_label} GWAS colocalizes with {args.gene} {args.cell_type.lower()} ieQTL", fontsize=8.4, fontweight="bold")
    fig.text(0.12, 0.912, rf"PP$_{{H4}}$ = {args.pph4:.3f}; PP$_{{H3}}$ = {args.pph3:.3f}", fontsize=7.2)
    fig.text(0.47, 0.912, f"Shared lead SNP: {args.lead}", fontsize=7.2)
    fig.text(0.12, 0.045, "Points are colored by LD with the shared lead SNP; diamond indicates the lead SNP.", fontsize=6.4, color="#555555")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = f"{args.trait_label.replace(' ', '_')}_{args.tissue_label.replace(' ', '_')}_{args.cell_type.replace(' ', '_')}_{args.gene}_GWAS_ieQTL_coloc"
    fig.savefig(outdir / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(outdir / f"{stem}.tiff", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    pd.DataFrame(
        [{
            "trait": args.trait_label,
            "tissue": args.tissue_label,
            "cell_type": args.cell_type,
            "gene": args.gene,
            "shared_lead_snp": args.lead,
            "ieqtl_lead_snp": ieqtl_lead["variant_id"],
            "ieqtl_lead_pval_g2": ieqtl_lead["pval_g2"],
            "PPH3": args.pph3,
            "PPH4": args.pph4,
            "n_shared_snps": args.shared_snps,
        }]
    ).to_csv(outdir / f"{stem}.source_data.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
