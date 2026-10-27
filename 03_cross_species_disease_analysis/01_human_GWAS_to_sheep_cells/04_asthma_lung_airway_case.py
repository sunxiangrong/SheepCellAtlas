#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import shutil
import tempfile
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse
from scipy.stats import ranksums, ttest_ind


ASTHMA_TRAIT = "categorical-20002-both_sexes-1111"
ASTHMA_TRAIT_CN = "哮喘"


def read_h5ad_gz(path: Path):
    if path.suffix == ".gz":
        tmpdir = tempfile.TemporaryDirectory()
        out = Path(tmpdir.name) / path.name.replace(".gz", "")
        with gzip.open(path, "rb") as f_in, open(out, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        adata = sc.read_h5ad(out)
        tmpdir.cleanup()
        return adata
    return sc.read_h5ad(path)


def bh_fdr(pvals: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvals, errors="coerce").fillna(1.0).clip(0, 1)
    n = len(p)
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return pd.Series(out, index=p.index)


def row_zscore(df: pd.DataFrame) -> pd.DataFrame:
    arr = df.to_numpy(dtype=float)
    mean = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, keepdims=True)
    keep = (sd[:, 0] > 1e-8) & np.isfinite(sd[:, 0])
    arr = (arr[keep] - mean[keep]) / sd[keep]
    return pd.DataFrame(arr, index=df.index[keep], columns=df.columns)


def normalize_counts(x):
    """Apply the established count-scale heuristic without densifying sparse data."""
    if sparse.issparse(x):
        x = x.tocsr().astype(np.float32, copy=False)
        values = x.data
        maximum = float(values.max()) if values.size else 0.0
        minimum = float(values.min()) if values.size else 0.0
        if minimum < 0:
            return x
        if maximum > 50:
            library_size = np.asarray(x.sum(axis=1)).ravel()
            scale = np.divide(
                10_000.0,
                library_size,
                out=np.zeros_like(library_size, dtype=np.float32),
                where=library_size > 0,
            )
            x = sparse.diags(scale).dot(x).tocsr()
            x.data = np.log1p(x.data)
            return x
        if maximum > 10:
            x.data = np.log1p(x.data)
        return x

    x = np.asarray(x, dtype=np.float32)
    if x.size == 0:
        return x
    if np.nanmax(x) > 50:
        library_size = x.sum(axis=1, keepdims=True)
        library_size[library_size == 0] = 1
        return np.log1p(x / library_size * 10_000)
    if np.nanmax(x) > 10:
        return np.log1p(x)
    return x


def pseudobulk_by(adata, group_cols: list[str], celltype: str | None = None) -> pd.DataFrame:
    obs = adata.obs.copy()
    if celltype is not None:
        keep = obs["cell_set_clean"].astype(str).eq(celltype).to_numpy()
        adata = adata[keep].copy()
        obs = adata.obs.copy()
    groups = obs[group_cols].astype(str).agg("__".join, axis=1)
    x = adata.X
    x = normalize_counts(x)
    rows = {}
    group_values = groups.to_numpy()
    for g in sorted(groups.unique()):
        idx = np.flatnonzero(group_values == g)
        if len(idx) >= 20:
            mean = x[idx].mean(axis=0)
            rows[g] = np.asarray(mean).ravel()
    return pd.DataFrame(rows, index=adata.var_names.astype(str))


def load_one2one(
    ortholog_path: Path,
    sheep_genes: set[str],
    human_genes: set[str],
    human_gene_column: str,
) -> tuple[pd.DataFrame, str]:
    orth = pd.read_csv(ortholog_path, sep="\t")
    required = {"Human_homology_type", "Sheep_gene_name", human_gene_column}
    missing = required.difference(orth.columns)
    if missing:
        raise ValueError(f"Ortholog table is missing columns: {sorted(missing)}")
    orth = orth[
        orth["Human_homology_type"].eq("ortholog_one2one")
        & orth["Sheep_gene_name"].notna()
        & orth[human_gene_column].notna()
    ].copy()
    common = orth[
        orth["Sheep_gene_name"].isin(sheep_genes)
        & orth[human_gene_column].isin(human_genes)
    ].copy()
    common = common.drop_duplicates("Sheep_gene_name").drop_duplicates(human_gene_column)
    return common, human_gene_column


def make_conserved_modules(
    sheep_pb: pd.DataFrame,
    human_cell_pb: pd.DataFrame,
    ortholog_path: Path,
    out_dir: Path,
    human_gene_column: str,
) -> pd.DataFrame:
    sheep = sheep_pb.drop_duplicates("gene").set_index("gene")
    human = human_cell_pb.copy()
    common, human_col = load_one2one(
        ortholog_path, set(sheep.index), set(human.index), human_gene_column
    )
    sheep_m = sheep.loc[common["Sheep_gene_name"].values]
    human_m = human.loc[common[human_col].values]
    sheep_m.index = common[human_col].values
    human_m.index = common[human_col].values
    shared = sheep_m.index.intersection(human_m.index)
    sheep_z = row_zscore(sheep_m.loc[shared])
    human_z = row_zscore(human_m.loc[shared])
    shared = sheep_z.index.intersection(human_z.index)
    sheep_z = sheep_z.loc[shared]
    human_z = human_z.loc[shared]

    epithelial_sheep = [c for c in ["Club cells", "Ciliated cells", "Alveolar epithelial type 1 cells", "Alveolar epithelial type 2 cells"] if c in sheep_z.columns]
    modules = pd.DataFrame(index=shared)
    modules["gene"] = shared
    modules["epithelial_conserved_module"] = False
    modules["cd4_conserved_module"] = False
    if "Epithelial cells" in human_z.columns and epithelial_sheep:
        modules["epithelial_conserved_module"] = (
            (human_z["Epithelial cells"] >= 0.75)
            & (sheep_z[epithelial_sheep].max(axis=1) >= 0.75)
        ).to_numpy()
    if "CD4 T cells" in human_z.columns and "CD4+ T cells" in sheep_z.columns:
        modules["cd4_conserved_module"] = (
            (human_z["CD4 T cells"] >= 0.75)
            & (sheep_z["CD4+ T cells"] >= 0.75)
        ).to_numpy()
    modules["asthma_relevant_conserved_module"] = modules["epithelial_conserved_module"] | modules["cd4_conserved_module"]
    modules.to_csv(out_dir / "asthma_lung_conserved_gene_modules.tsv", sep="\t", index=False)
    return modules


def plot_scpagwas_dot(human_tissue_tsv: Path, out_dir: Path) -> pd.DataFrame:
    df = pd.read_csv(human_tissue_tsv, sep="\t")
    d = df[(df["trait"].eq(ASTHMA_TRAIT)) & (df["context"].eq("Lung"))].copy()
    d["mean_TRS"] = pd.to_numeric(d["mean_TRS"], errors="coerce")
    d["fdr_or_adjp"] = pd.to_numeric(d["fdr_or_adjp"], errors="coerce")
    d["neglog10_fdr"] = -np.log10(d["fdr_or_adjp"].clip(lower=1e-300))
    d = d.sort_values(["mean_TRS", "neglog10_fdr"], ascending=[False, False])
    d.to_csv(out_dir / "asthma_lung_scpagwas_TRS_source.tsv", sep="\t", index=False)

    sns.set_theme(style="white", context="paper", font_scale=0.95)
    fig, ax = plt.subplots(figsize=(6.8, 2.7))
    sizes = 45 + d["mean_TRS"].clip(lower=0).fillna(0) * 190
    sc = ax.scatter(
        np.arange(len(d)),
        d["mean_TRS"],
        s=sizes,
        c=d["neglog10_fdr"],
        cmap="Reds",
        edgecolor="#333333",
        linewidth=0.35,
        alpha=0.9,
    )
    ax.axhline(0.5, color="#555555", ls="--", lw=0.8)
    ax.set_xticks(np.arange(len(d)))
    # Use English labels here because the IAS plotting environment lacks a CJK font.
    ax.set_xticklabels(d["cell_group"], rotation=35, ha="right")
    ax.set_ylabel("mean TRS")
    ax.set_xlabel("")
    ax.set_ylim(min(0, d["mean_TRS"].min() - 0.08), d["mean_TRS"].max() + 0.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    cbar = fig.colorbar(sc, ax=ax, pad=0.015)
    cbar.set_label("-log10(FDR)")
    fig.savefig(out_dir / "panelB_asthma_lung_scpagwas_trs_dot.png", dpi=450, bbox_inches="tight")
    fig.savefig(out_dir / "panelB_asthma_lung_scpagwas_trs_dot.pdf", bbox_inches="tight")
    plt.close(fig)
    return d


def plot_correlation_heatmap(corr_tsv: Path, out_dir: Path) -> None:
    corr = pd.read_csv(corr_tsv, sep="\t", index_col=0)
    sns.set_theme(style="white", context="paper", font_scale=0.9)
    cg = sns.clustermap(
        corr.astype(float),
        cmap="vlag",
        center=0,
        linewidths=0.25,
        linecolor="#F0F0F0",
        figsize=(7.2, 6.5),
        cbar_kws={"label": "Spearman correlation"},
        xticklabels=True,
        yticklabels=True,
    )
    cg.ax_heatmap.set_xlabel("Human airway cell type")
    cg.ax_heatmap.set_ylabel("Sheep lung cell type")
    cg.fig.savefig(out_dir / "panelC_sheep_lung_human_airway_celltype_correlation.png", dpi=450, bbox_inches="tight")
    cg.fig.savefig(out_dir / "panelC_sheep_lung_human_airway_celltype_correlation.pdf", bbox_inches="tight")
    plt.close(cg.fig)


def deg_for_celltype(adata, celltype: str) -> pd.DataFrame:
    pb = pseudobulk_by(adata, ["phenotype", "sample"], celltype=celltype)
    aa_cols = [c for c in pb.columns if c.startswith("AA__")]
    ana_cols = [c for c in pb.columns if c.startswith("ANA__")]
    aa = pb[aa_cols]
    ana = pb[ana_cols]
    logfc = aa.mean(axis=1) - ana.mean(axis=1)
    stat, pvals = ttest_ind(aa.T, ana.T, axis=0, equal_var=False, nan_policy="omit")
    out = pd.DataFrame(
        {
            "gene": pb.index,
            "celltype": celltype,
            "logFC_AA_vs_ANA": logfc.values,
            "abs_logFC": np.abs(logfc.values),
            "pvalue": np.nan_to_num(pvals, nan=1.0),
            "n_AA_pseudobulk": len(aa_cols),
            "n_ANA_pseudobulk": len(ana_cols),
        }
    )
    out["FDR"] = bh_fdr(out["pvalue"])
    return out.sort_values(["FDR", "abs_logFC"], ascending=[True, False])


def celllevel_response_for_celltype(
    adata,
    celltype: str,
    case_group: str = "AA",
    stimulated: str = "Ag",
    baseline: str = "Pre",
) -> pd.DataFrame:
    obs = adata.obs.copy()
    keep = (
        obs["cell_set_clean"].astype(str).eq(celltype)
        & obs["phenotype"].astype(str).eq(case_group)
        & obs["sample"].astype(str).isin([stimulated, baseline])
    ).to_numpy()
    sub = adata[keep].copy()
    group = sub.obs["sample"].astype(str).to_numpy()
    # Keep the full cell-by-gene matrix sparse. Rank tests need only one gene
    # at a time, so densify two short vectors rather than the entire matrix.
    x = normalize_counts(sub.X)
    stim_idx = np.where(group == stimulated)[0]
    base_idx = np.where(group == baseline)[0]
    stim = x[stim_idx, :]
    base = x[base_idx, :]
    logfc = np.asarray(stim.mean(axis=0)).ravel() - np.asarray(base.mean(axis=0)).ravel()
    pvals = np.ones(x.shape[1], dtype=float)
    for j in range(x.shape[1]):
        try:
            stim_values = stim[:, j].toarray().ravel() if sparse.issparse(stim) else stim[:, j]
            base_values = base[:, j].toarray().ravel() if sparse.issparse(base) else base[:, j]
            pvals[j] = ranksums(stim_values, base_values).pvalue
        except Exception:
            pvals[j] = 1.0
    out = pd.DataFrame(
        {
            "gene": sub.var_names.astype(str),
            "celltype": celltype,
            "contrast": f"{case_group}_{stimulated}_vs_{baseline}",
            "logFC_Ag_vs_Pre_in_AA": logfc,
            "abs_logFC": np.abs(logfc),
            "pvalue": np.nan_to_num(pvals, nan=1.0),
            "n_Ag_cells": len(stim_idx),
            "n_Pre_cells": len(base_idx),
        }
    )
    out["FDR"] = bh_fdr(out["pvalue"])
    return out.sort_values(["FDR", "abs_logFC"], ascending=[True, False])


def plot_volcano(deg: pd.DataFrame, modules: pd.DataFrame, module_col: str, label: str, out_dir: Path) -> pd.DataFrame:
    d = deg.merge(modules[["gene", module_col]], on="gene", how="left")
    d[module_col] = d[module_col].fillna(False).astype(bool)
    d["class"] = "Other genes"
    fc_col = "logFC_AA_vs_ANA" if "logFC_AA_vs_ANA" in d.columns else "logFC_Ag_vs_Pre_in_AA"
    contrast_label = "asthma / control" if fc_col == "logFC_AA_vs_ANA" else "allergen / baseline in asthma"
    sig = (d["FDR"] <= 0.05) & (d["abs_logFC"] >= 0.25)
    d.loc[sig & d[module_col], "class"] = "Conserved-module DEG"
    d.loc[sig & ~d[module_col], "class"] = "Other DEG"
    d["neglog10p"] = -np.log10(d["pvalue"].clip(lower=1e-300))
    source_prefix = "panelD_asthma_Ag_vs_Pre_in_AA" if "logFC_Ag_vs_Pre_in_AA" in d.columns else "panelD_asthma_AA_vs_ANA"
    d.to_csv(out_dir / f"{source_prefix}_{label}.volcano_source.tsv", sep="\t", index=False)

    palette = {"Other genes": "#D2D2D2", "Other DEG": "#E6862F", "Conserved-module DEG": "#3B6FA3"}
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    for cls in ["Other genes", "Other DEG", "Conserved-module DEG"]:
        sub = d[d["class"].eq(cls)]
        ax.scatter(
            sub[fc_col],
            sub["neglog10p"],
            s=7 if cls == "Other genes" else 18,
            c=palette[cls],
            alpha=0.25 if cls == "Other genes" else 0.85,
            linewidths=0,
            label=cls,
        )
    ax.axvline(0.25, ls="--", lw=0.7, color="#666666")
    ax.axvline(-0.25, ls="--", lw=0.7, color="#666666")
    ax.axhline(-np.log10(0.05), ls=":", lw=0.7, color="#666666")
    ymax = min(60, max(8, np.nanpercentile(d["neglog10p"], 99.0) * 1.08))
    ax.set_ylim(0, ymax)
    ax.set_xlabel(f"logFC ({contrast_label})")
    ax.set_ylabel("-log10(P)")
    ax.set_title(label.replace("_", " "), fontsize=10, fontweight="bold")
    ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    prefix = "panelD_asthma_Ag_vs_Pre_in_AA" if fc_col == "logFC_Ag_vs_Pre_in_AA" else "panelD_asthma_AA_vs_ANA"
    fig.savefig(out_dir / f"{prefix}_{label}.volcano.png", dpi=450, bbox_inches="tight")
    fig.savefig(out_dir / f"{prefix}_{label}.volcano.pdf", bbox_inches="tight")
    plt.close(fig)
    return d


def plot_conserved_fraction(volcano_tables: list[pd.DataFrame], out_dir: Path) -> None:
    rows = []
    for d in volcano_tables:
        celltype = d["celltype"].iloc[0]
        sig = d["class"].ne("Other genes")
        rows.append(
            {
                "celltype": celltype,
                "conserved": int((sig & d["class"].eq("Conserved-module DEG")).sum()),
                "other_deg": int((sig & d["class"].eq("Other DEG")).sum()),
            }
        )
    s = pd.DataFrame(rows)
    s["total"] = s["conserved"] + s["other_deg"]
    s["fraction_conserved"] = s["conserved"] / s["total"].replace(0, np.nan)
    s.to_csv(out_dir / "panelE_asthma_conserved_vs_other_deg_summary.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(4.2, 2.2))
    y = np.arange(len(s))
    ax.barh(y, s["conserved"], color="#3B6FA3", label="Conserved-module DEG")
    ax.barh(y, s["other_deg"], left=s["conserved"], color="#E6862F", label="Other DEG")
    ax.set_yticks(y)
    ax.set_yticklabels(s["celltype"])
    ax.set_xlabel("Number of DEG")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.savefig(out_dir / "panelE_asthma_conserved_vs_other_deg_bar.png", dpi=450, bbox_inches="tight")
    fig.savefig(out_dir / "panelE_asthma_conserved_vs_other_deg_bar.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--human-tissue-tsv", required=True)
    parser.add_argument("--corr-tsv", required=True)
    parser.add_argument("--sheep-pseudobulk", required=True)
    parser.add_argument("--human-h5ad", required=True)
    parser.add_argument("--ortholog", required=True)
    parser.add_argument(
        "--human-ortholog-column",
        choices=["Human_gene_name_hg19", "Human_gene_name_hg38"],
        default="Human_gene_name_hg19",
        help="Human gene-symbol column matching the supplied human expression object.",
    )
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["axes.unicode_minus"] = False

    plot_scpagwas_dot(Path(args.human_tissue_tsv), out_dir)
    plot_correlation_heatmap(Path(args.corr_tsv), out_dir)

    adata = read_h5ad_gz(Path(args.human_h5ad))
    human_cell_pb = pseudobulk_by(adata, ["cell_set_clean"])
    sheep_pb = pd.read_csv(args.sheep_pseudobulk, sep="\t")
    modules = make_conserved_modules(
        sheep_pb,
        human_cell_pb,
        Path(args.ortholog),
        out_dir,
        args.human_ortholog_column,
    )

    epithelial_deg = celllevel_response_for_celltype(adata, "Epithelial cells")
    cd4_deg = celllevel_response_for_celltype(adata, "CD4 T cells")
    epithelial_table = plot_volcano(epithelial_deg, modules, "epithelial_conserved_module", "Epithelial_cells", out_dir)
    cd4_table = plot_volcano(cd4_deg, modules, "cd4_conserved_module", "CD4_T_cells", out_dir)
    plot_conserved_fraction([epithelial_table, cd4_table], out_dir)
    print("WROTE", out_dir)


if __name__ == "__main__":
    main()
