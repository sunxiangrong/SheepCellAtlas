from __future__ import annotations

import argparse
import gzip
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt


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


def pseudobulk_from_h5ad(path: Path, label: str, layer: str | None = None, query: str | None = None) -> pd.Series:
    adata = read_h5ad_gz(path)
    if query:
        keep = adata.obs.eval(query)
        adata = adata[keep.to_numpy()].copy()
    if layer:
        if layer == "raw_adata":
            x = adata.raw.X
            genes = adata.raw.var_names.astype(str)
        else:
            x = adata.layers[layer]
            genes = adata.var_names.astype(str)
    else:
        x = adata.X
        genes = adata.var_names.astype(str)
    if not isinstance(x, np.ndarray):
        avg = np.asarray(x.mean(axis=0)).ravel()
    else:
        avg = x.mean(axis=0)
    return pd.Series(avg, index=genes, name=label)


def row_zscore(df: pd.DataFrame) -> pd.DataFrame:
    arr = df.to_numpy(dtype=float)
    mean = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, keepdims=True)
    keep = (sd[:, 0] > 1e-8) & np.isfinite(sd[:, 0])
    arr = (arr[keep] - mean[keep]) / sd[keep]
    return pd.DataFrame(arr, index=df.index[keep], columns=df.columns)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sheep-tissue-pseudobulk", required=True)
    parser.add_argument("--human-h5ad", action="append", required=True, help="label:path:layer or label:path")
    parser.add_argument("--ortholog", required=True)
    parser.add_argument(
        "--human-ortholog-column",
        choices=["Human_gene_name_hg19", "Human_gene_name_hg38"],
        default="Human_gene_name_hg19",
        help="Human gene-symbol column that matches the supplied H5AD objects.",
    )
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--title", default="Sheep-human tissue pseudobulk correlation")
    args = parser.parse_args()

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    sheep = pd.read_csv(args.sheep_tissue_pseudobulk, sep="\t").drop_duplicates("gene").set_index("gene")

    human_series = []
    for item in args.human_h5ad:
        parts = item.split(":", 2)
        if len(parts) < 2:
            raise ValueError("--human-h5ad must be label:path[:layer]")
        label, path = parts[0], parts[1]
        layer = parts[2] if len(parts) == 3 and parts[2] else None
        human_series.append(pseudobulk_from_h5ad(Path(path), label, layer=layer))
    human = pd.concat(human_series, axis=1)

    orth = pd.read_csv(args.ortholog, sep="\t")
    human_gene_col = args.human_ortholog_column
    required_ortholog_columns = {
        "Human_homology_type", "Sheep_gene_name", human_gene_col
    }
    missing = required_ortholog_columns.difference(orth.columns)
    if missing:
        raise ValueError(f"Ortholog table is missing columns: {sorted(missing)}")
    orth = orth[
        (orth["Human_homology_type"] == "ortholog_one2one")
        & orth["Sheep_gene_name"].notna()
        & orth[human_gene_col].notna()
    ].drop_duplicates(["Sheep_gene_name", human_gene_col])
    common = orth[
        orth["Sheep_gene_name"].isin(sheep.index)
        & orth[human_gene_col].isin(human.index)
    ].copy()
    common = common.drop_duplicates("Sheep_gene_name").drop_duplicates(human_gene_col)

    sheep_m = sheep.loc[common["Sheep_gene_name"].values]
    human_m = human.loc[common[human_gene_col].values]
    sheep_m.index = common[human_gene_col].values
    human_m.index = common[human_gene_col].values

    valid = np.isfinite(sheep_m.to_numpy()).all(axis=1) & np.isfinite(human_m.to_numpy()).all(axis=1)
    sheep_m = sheep_m.loc[valid]
    human_m = human_m.loc[valid]
    sheep_z = row_zscore(sheep_m)
    human_z = row_zscore(human_m)
    shared = sheep_z.index.intersection(human_z.index)
    sheep_z = sheep_z.loc[shared]
    human_z = human_z.loc[shared]

    corr = pd.DataFrame(index=sheep_z.columns, columns=human_z.columns, dtype=float)
    for s in corr.index:
        for h in corr.columns:
            corr.loc[s, h] = pd.Series(sheep_z[s].values).corr(pd.Series(human_z[h].values), method="spearman")
    corr = corr.replace([np.inf, -np.inf], np.nan).fillna(0)
    corr.to_csv(str(out_prefix) + ".correlation.tsv", sep="\t")
    common.to_csv(str(out_prefix) + ".ortholog_used.tsv", sep="\t", index=False)

    sns.set_theme(style="white", font_scale=0.8)
    width = max(7, 0.38 * corr.shape[1] + 4)
    height = max(8, 0.25 * corr.shape[0] + 3)
    cg = sns.clustermap(
        corr.astype(float),
        cmap="vlag",
        center=0,
        linewidths=0.15,
        linecolor="#eeeeee",
        figsize=(width, height),
        cbar_kws={"label": "Spearman correlation"},
        xticklabels=True,
        yticklabels=True,
    )
    cg.fig.suptitle(f"{args.title}\n{len(shared)} one-to-one ortholog genes", y=1.02, fontsize=11)
    cg.ax_heatmap.set_xlabel("Human tissue/dataset pseudobulk")
    cg.ax_heatmap.set_ylabel("Sheep tissue pseudobulk")
    cg.fig.savefig(str(out_prefix) + ".heatmap.png", dpi=300, bbox_inches="tight")
    cg.fig.savefig(str(out_prefix) + ".heatmap.pdf", bbox_inches="tight")
    print("WROTE", str(out_prefix) + ".heatmap.png")
    print("ortholog_genes", len(shared), "sheep_tissues", corr.shape[0], "human_tissues", corr.shape[1])


if __name__ == "__main__":
    main()
