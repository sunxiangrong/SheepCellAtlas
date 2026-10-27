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
from scipy import sparse


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


def human_pseudobulk(adata, group_col: str, min_cells: int = 20, layer: str | None = None) -> pd.DataFrame:
    groups = adata.obs[group_col].astype(str)
    keep_groups = groups.value_counts()
    keep_groups = keep_groups[keep_groups >= min_cells].index.tolist()
    use = groups.isin(keep_groups).to_numpy()
    adata = adata[use].copy()
    groups = adata.obs[group_col].astype(str)
    if layer:
        if layer == "raw_adata":
            if adata.raw is None:
                raise ValueError("Requested raw_adata but adata.raw is None")
            x = adata.raw.X
            var_names = adata.raw.var_names.astype(str)
        else:
            if layer not in adata.layers:
                raise ValueError(f"Layer not found: {layer}; available={list(adata.layers.keys())}")
            x = adata.layers[layer]
            var_names = adata.var_names.astype(str)
    else:
        x = adata.X
        var_names = adata.var_names.astype(str)
    if sparse.issparse(x):
        # Keep sparse matrices sparse: human single-cell objects can be too
        # large to densify safely in memory.
        x = x.tocsr().astype(np.float32, copy=False)
        values = x.data
        xmin = float(values.min()) if values.size else 0.0
        xmax = float(values.max()) if values.size else 0.0
        if xmin >= 0 and xmax > 50:
            library_size = np.asarray(x.sum(axis=1)).ravel()
            scale = np.divide(
                10_000.0,
                library_size,
                out=np.zeros_like(library_size, dtype=np.float32),
                where=library_size > 0,
            )
            x = sparse.diags(scale).dot(x).tocsr()
            x.data = np.log1p(x.data)
        elif xmin >= 0 and xmax > 10:
            x.data = np.log1p(x.data)
    else:
        x = np.asarray(x, dtype=np.float32)
        xmin = np.nanmin(x)
        xmax = np.nanmax(x)
        if xmin >= 0 and xmax > 50:
            library_size = x.sum(axis=1, keepdims=True)
            library_size[library_size == 0] = 1
            x = np.log1p(x / library_size * 10_000)
        elif xmin >= 0 and xmax > 10:
            x = np.log1p(x)
    rows = {}
    group_values = groups.to_numpy()
    for g in sorted(set(groups)):
        idx = np.flatnonzero(group_values == g)
        rows[g] = np.asarray(x[idx].mean(axis=0)).ravel()
    pb = pd.DataFrame(rows, index=var_names)
    pb.index.name = "human_gene"
    return pb


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sheep-pseudobulk", required=True)
    parser.add_argument("--human-h5ad", required=True)
    parser.add_argument("--human-group-col", required=True)
    parser.add_argument("--human-layer", default="")
    parser.add_argument("--ortholog", required=True)
    parser.add_argument(
        "--human-ortholog-column",
        choices=["Human_gene_name_hg19", "Human_gene_name_hg38"],
        default="Human_gene_name_hg19",
        help="Human gene-symbol column that matches the supplied H5AD object.",
    )
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--title", default="Sheep-human cell type correlation")
    args = parser.parse_args()

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    sheep = pd.read_csv(args.sheep_pseudobulk, sep="\t")
    sheep = sheep.drop_duplicates("gene").set_index("gene")

    adata = read_h5ad_gz(Path(args.human_h5ad))
    human = human_pseudobulk(adata, args.human_group_col, layer=args.human_layer or None)

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
    ].copy()
    orth = orth.drop_duplicates(["Sheep_gene_name", human_gene_col])

    common = orth[
        orth["Sheep_gene_name"].isin(sheep.index)
        & orth[human_gene_col].isin(human.index)
    ].copy()

    common = common.drop_duplicates(human_gene_col).drop_duplicates("Sheep_gene_name")
    sheep_m = sheep.loc[common["Sheep_gene_name"].values]
    human_m = human.loc[common[human_gene_col].values]
    sheep_m.index = common[human_gene_col].values
    human_m.index = common[human_gene_col].values
    valid_genes = np.isfinite(sheep_m.to_numpy()).all(axis=1) & np.isfinite(human_m.to_numpy()).all(axis=1)
    sheep_m = sheep_m.loc[valid_genes]
    human_m = human_m.loc[valid_genes]
    # Compare cell-type-specific patterns rather than global housekeeping expression:
    # z-score each gene across cell types within each species, then correlate profiles.
    def row_zscore(df: pd.DataFrame) -> pd.DataFrame:
        arr = df.to_numpy(dtype=float)
        mean = np.nanmean(arr, axis=1, keepdims=True)
        sd = np.nanstd(arr, axis=1, keepdims=True)
        keep = (sd[:, 0] > 1e-8) & np.isfinite(sd[:, 0])
        arr = (arr[keep] - mean[keep]) / sd[keep]
        return pd.DataFrame(arr, index=df.index[keep], columns=df.columns)

    sheep_m = row_zscore(sheep_m)
    human_m = row_zscore(human_m)
    shared_idx = sheep_m.index.intersection(human_m.index)
    sheep_m = sheep_m.loc[shared_idx]
    human_m = human_m.loc[shared_idx]

    corr = pd.DataFrame(index=sheep_m.columns, columns=human_m.columns, dtype=float)
    for s in sheep_m.columns:
        for h in human_m.columns:
            corr.loc[s, h] = pd.Series(sheep_m[s].values).corr(pd.Series(human_m[h].values), method="spearman")
    corr = corr.replace([np.inf, -np.inf], np.nan)
    corr = corr.dropna(axis=0, how="all").dropna(axis=1, how="all")
    corr = corr.fillna(0)

    corr.to_csv(str(out_prefix) + ".correlation.tsv", sep="\t")
    common.to_csv(str(out_prefix) + ".ortholog_used.tsv", sep="\t", index=False)

    n_genes = len(common)
    width = max(8, 0.45 * corr.shape[1] + 3)
    height = max(7, 0.35 * corr.shape[0] + 2)
    sns.set_theme(style="white", font_scale=0.75)
    cg = sns.clustermap(
        corr.astype(float),
        cmap="vlag",
        center=0,
        linewidths=0.2,
        linecolor="#eeeeee",
        figsize=(width, height),
        cbar_kws={"label": "Spearman correlation"},
        xticklabels=True,
        yticklabels=True,
    )
    cg.fig.suptitle(f"{args.title}\n{n_genes} one-to-one ortholog genes", y=1.02, fontsize=11)
    cg.ax_heatmap.set_xlabel("Human cell type")
    cg.ax_heatmap.set_ylabel("Sheep cell type")
    cg.fig.savefig(str(out_prefix) + ".heatmap.png", dpi=300, bbox_inches="tight")
    cg.fig.savefig(str(out_prefix) + ".heatmap.pdf", bbox_inches="tight")
    plt.close(cg.fig)

    print("WROTE", str(out_prefix) + ".correlation.tsv")
    print("WROTE", str(out_prefix) + ".heatmap.png")
    print("ortholog_genes", n_genes, "sheep_groups", corr.shape[0], "human_groups", corr.shape[1])


if __name__ == "__main__":
    main()
