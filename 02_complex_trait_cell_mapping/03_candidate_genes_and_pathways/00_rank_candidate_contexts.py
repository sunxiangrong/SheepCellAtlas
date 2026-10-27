#!/usr/bin/env python3
"""Plot ranked PCC curves for ten pre-specified sheep scPagwas contexts."""

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection


CONTEXTS = [
    ("Spleen", "ALLGE3", "B cells", "fav_tissue_ds2000_Spleen_ALLGE3_B_cells.tsv"),
    ("Uterus", "ALLGE3", "gamma-delta T cells", "fav_tissue_ds2000_Uterus_ALLGE3_T_cells.tsv"),
    ("Oviduct", "ALLGE3", "Luminal cells", "fav_tissue_ds2000_Oviduct_ALLGE3_Luminal_cells.tsv"),
    ("Conceptus", "ALLGE3", "Hypoblasts cells", "fav_tissue_ds2000_Conceptus_ALLGE3_Hypoblasts_cells.tsv"),
    ("Lymph node", "YELN", "T cells", "fav_tissue_ds2000_Lymphnode_YELN_T_cells.tsv"),
    ("Ovary", "CH", "NK T cells", "fav_tissue_ds2000_Ovary_CH_NK_T_cells.tsv"),
    ("Thymus", "CH", "T cells", "fav_tissue_ds2000_Thymus_CH_T_cells.tsv"),
    ("Lung", "CD", "CD4+ T cells", "fav_tissue_ds2000_Lung_CD_CD4_T_cells.tsv"),
    ("Rectum", "EL", "C7+ stromal cells", "fav_tissue_ds2000_Rectum_EL_C7_stromal_cells_C7_SCs.tsv"),
    ("Reticulum", "NL", "Spinous cells", "fav_tissue_ds2000_Reticulum_NL_Spinous_cells.tsv"),
]


def keep_gene(gene: str) -> bool:
    s = str(gene).strip().upper()
    if s.startswith(("LOC", "RPS", "RPL", "KEF53")):
        return False
    if s.startswith(("MT-", "MT.", "MT_", "MTND", "MTCO", "MTATP", "MTCYB", "MTTL", "MTRNR", "MITO")):
        return False
    return s not in {"ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "CYTB", "ATP6", "ATP8", "RNR1", "RNR2"}


def load_context(input_dir: Path, filename: str, include_all_genes=False) -> pd.DataFrame:
    path = input_dir / filename
    df = pd.read_csv(path, sep="\t")
    if not {"gene", "PCC_clean"}.issubset(df.columns):
        raise ValueError(f"Missing gene/PCC_clean columns: {path}")
    if not include_all_genes:
        df = df[df["gene"].map(keep_gene)].copy()
    df["PCC_clean"] = pd.to_numeric(df["PCC_clean"], errors="coerce")
    df = df.dropna(subset=["PCC_clean"]).sort_values("PCC_clean", ascending=False).reset_index(drop=True)
    df["order"] = np.arange(1, len(df) + 1)
    if df.empty or not np.isfinite(df["PCC_clean"]).all():
        raise ValueError(f"Empty or non-finite PCC values: {path}")
    return df


def add_labels(ax, df: pd.DataFrame) -> None:
    top = df.nlargest(3, "PCC_clean").copy()
    x_max = max(1, len(df))
    y_span = max(0.08, float(df["PCC_clean"].max() - df["PCC_clean"].min()))
    # Place labels close to the right tail while staggering their vertical positions.
    for i, (_, row) in enumerate(top.iterrows()):
        x = float(row["order"])
        y = float(row["PCC_clean"])
        # The x-axis is reversed to match the reference: rank 1 is at the right.
        tx = min(x_max, x + 0.17 * x_max)
        ty = y + (0.08 - i * 0.055) * y_span
        label_pcc = f"{y:.3f}" if abs(y) < 0.01 else f"{y:.2f}"
        ax.annotate(
            f"{row['gene']} (PCC = {label_pcc})",
            xy=(x, y),
            xytext=(tx, ty),
            textcoords="data",
            ha="right",
            va="bottom",
            fontsize=5.6,
            color="#A94738",
            arrowprops=dict(arrowstyle="-", lw=0.55, color="#A94738", shrinkA=1, shrinkB=2),
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True,
                        help="Directory containing pcc_top1000_<context>.tsv files")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--contexts", type=Path,
                        help="TSV columns: unit, trait, celltype, filename (relative to input-dir)")
    parser.add_argument("--include-all-genes", action="store_true",
                        help="Disable manuscript label-based LOC/ribosomal/mitochondrial filtering")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "font.size": 7,
    })
    cmap = LinearSegmentedColormap.from_list("pcc", ["#385A9B", "#EEEAF4", "#E76F51"], N=256)
    contexts = CONTEXTS
    if args.contexts:
        table = pd.read_csv(args.contexts, sep="\t", dtype=str)
        columns = ["unit", "trait", "celltype", "filename"]
        if not set(columns).issubset(table.columns) or table.empty or table[columns].isna().any().any():
            parser.error("Context TSV requires nonempty unit, trait, celltype, filename columns")
        contexts = list(table[columns].itertuples(index=False, name=None))
    else:
        contexts = [(u, t, c, "pcc_top1000_" + f) for u, t, c, f in CONTEXTS]
    frames = [load_context(args.input_dir, filename, args.include_all_genes) for _, _, _, filename in contexts]
    # Retain a reusable, auditable result table independently of figure output.
    pd.concat([df.assign(unit=u, trait=t, celltype=c, source_file=f)
               for (u, t, c, f), df in zip(contexts, frames)], ignore_index=True).to_csv(
        args.out_dir / "ranked_candidates.tsv", sep="\t", index=False)
    global_vmax = max(0.05, max(float(np.nanmax(np.abs(df["PCC_clean"]))) for df in frames))
    global_norm = TwoSlopeNorm(vmin=-global_vmax, vcenter=0.0, vmax=global_vmax)
    nrows = (len(contexts) + 1) // 2
    fig, axes = plt.subplots(nrows, 2, figsize=(10.2, max(3, 2.4 * nrows)), squeeze=False)
    axes = axes.ravel()
    for unused in axes[len(contexts):]:
        unused.set_visible(False)

    for ax, (unit, trait, celltype, filename), df in zip(axes, contexts, frames):
        vals = df["PCC_clean"].to_numpy(float)
        colors = cmap(global_norm(vals))
        points = np.column_stack([df["order"].to_numpy(float), vals])
        segments = np.stack([points[:-1], points[1:]], axis=1)
        ax.add_collection(LineCollection(segments, colors=cmap(global_norm((vals[:-1] + vals[1:]) / 2)), linewidths=0.85, alpha=0.6, zorder=1))
        ax.scatter(df["order"], vals, c=colors, s=7, linewidths=0, alpha=0.95, zorder=2)
        ax.axhline(0, color="#B8B8C0", lw=0.45, zorder=0)
        add_labels(ax, df)
        ax.set_title(f"{unit} | {trait} | {celltype}", fontsize=8, fontweight="bold", pad=5)
        ax.set_xlabel("Genes ranked by PCC", fontsize=6.8)
        ax.set_ylabel("PCC", fontsize=7.2)
        ax.tick_params(labelsize=6.2, length=2.5, pad=2)
        ax.set_xlim(0.5, len(df) + 0.5)
        ax.invert_xaxis()
        ax.grid(False)

    sm = ScalarMappable(norm=global_norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.tolist(), fraction=0.018, pad=0.018, aspect=32)
    cbar.set_label("PCC", fontsize=7)
    cbar.ax.tick_params(labelsize=6, width=0.5, length=2)
    fig.suptitle("Ranked PCC profiles across selected sheep scPagwas contexts", fontsize=10, fontweight="bold", y=0.995)
    note = "All supplied gene labels retained." if args.include_all_genes else "Genes are ranked after removing LOC, RPS/RPL, KEF53 and explicitly mitochondrial gene labels."
    fig.text(0.5, 0.006, note, ha="center", fontsize=6.2, color="#666666")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.965, bottom=0.035, hspace=0.55, wspace=0.30)
    base = args.out_dir / "pcc_ranked_contexts"
    fig.savefig(base.with_suffix(".png"), dpi=450, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    print(base.with_suffix(".png"))
    print(base.with_suffix(".pdf"))
    print(base.with_suffix(".svg"))


if __name__ == "__main__":
    main()
