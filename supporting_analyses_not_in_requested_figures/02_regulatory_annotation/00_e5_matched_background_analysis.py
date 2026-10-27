#!/usr/bin/env python3
"""Compare PPH3/PPH4-high ieQTL lead SNPs with matched cis-testing background SNPs."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm


mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7.5, "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.75,
})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--leads", required=True)
    parser.add_argument("--ieqtl-root", required=True)
    parser.add_argument("--e5-bed", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260718)
    args = parser.parse_args()
    if args.permutations < 1:
        parser.error("--permutations must be at least 1")
    return args


def wilson_interval(successes: int, total: int) -> tuple[float, float]:
    z = norm.ppf(0.975)
    p = successes / total
    denominator = 1 + z * z / total
    centre = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * total)) / total) / denominator
    return centre - half, centre + half


def load_intervals(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    data = pd.read_csv(path, sep="\t", compression="infer", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
    return {
        str(chrom): (
            block.sort_values("start").start.to_numpy(dtype=np.int64),
            block.sort_values("start").end.to_numpy(dtype=np.int64),
        )
        for chrom, block in data.groupby("chrom", sort=False)
    }


def overlap(frame: pd.DataFrame, intervals: dict[str, tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
    # Permutation draws may contain the same nominal-SNP row more than once.
    # Use positional row indices so repeated source indices are counted separately.
    frame = frame.reset_index(drop=True)
    result = np.zeros(len(frame), dtype=bool)
    for chrom, indexes in frame.groupby("chrom", sort=False).groups.items():
        if chrom not in intervals:
            continue
        starts, ends = intervals[chrom]
        positions = frame.loc[indexes, "pos"].to_numpy(dtype=np.int64) - 1
        at = np.searchsorted(starts, positions, side="right") - 1
        valid = at >= 0
        local = np.zeros(len(positions), dtype=bool)
        local[valid] = positions[valid] < ends[at[valid]]
        result[frame.index.get_indexer(indexes)] = local
    return result


def attach_nominal_and_background(leads: pd.DataFrame, root: Path) -> tuple[pd.DataFrame, dict[int, pd.DataFrame]]:
    annotated: list[pd.DataFrame] = []
    background_choices: dict[int, pd.DataFrame] = {}
    sequence = 0
    for (tissue, cell_type), block in leads.groupby(["tissue", "cell_type"], sort=False):
        requested = set(block["phenotype2"].astype(str))
        path = root / str(tissue) / f"{cell_type}.cis_qtl.tsv.gz"
        values: list[pd.DataFrame] = []
        for chunk in pd.read_csv(path, sep="\t", compression="gzip",
                                 usecols=["pheno_id", "variant_id", "start_distance", "af", "pval_g2"], chunksize=250_000):
            hit = chunk[chunk.pheno_id.astype(str).isin(requested)]
            if not hit.empty:
                values.append(hit)
        if not values:
            continue
        nominal = pd.concat(values, ignore_index=True)
        for lead in block.itertuples(index=False):
            gene_data = nominal[nominal.pheno_id.astype(str) == str(lead.phenotype2)].copy()
            actual = gene_data[gene_data.variant_id == lead.variant_id]
            if actual.empty:
                continue
            actual = actual.iloc[0]
            candidates = gene_data[(gene_data.pval_g2 > 0.05) & (gene_data.variant_id != lead.variant_id)].copy()
            if candidates.empty:
                continue
            lead_af = float(actual.af)
            lead_dist = math.log1p(abs(float(actual.start_distance)))
            candidates["match_score"] = (
                np.abs(candidates.af.astype(float) - lead_af) / 0.03
                + np.abs(np.log1p(np.abs(candidates.start_distance.astype(float))) - lead_dist) / 0.50
            )
            pool = candidates.nsmallest(min(100, len(candidates)), "match_score")
            record = lead._asdict()
            record.update({"af": actual.af, "start_distance": actual.start_distance, "lead_index": sequence})
            annotated.append(pd.DataFrame([record]))
            background_choices[sequence] = pool.reset_index(drop=True)
            sequence += 1
    if not annotated:
        raise ValueError("No lead variants matched the nominal ieQTL inputs and background criteria")
    return pd.concat(annotated, ignore_index=True), background_choices


def build_background_distribution(leads: pd.DataFrame, choices: dict[int, pd.DataFrame], intervals: dict[str, tuple[np.ndarray, np.ndarray]], permutations: int, seed: int) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    distributions: dict[str, list[float]] = {group: [] for group in ["PPH3-high", "PPH4-high"]}
    for _ in range(permutations):
        for group in distributions:
            selected = leads[leads.group == group]
            draws = []
            for row in selected.itertuples(index=False):
                pool = choices[row.lead_index]
                draws.append(pool.iloc[rng.integers(len(pool))])
            draw = pd.DataFrame(draws)
            chrom_pos = draw.variant_id.astype(str).str.split("_", n=1, expand=True)
            draw["chrom"] = "chr" + chrom_pos[0].str.replace("chr", "", regex=False)
            draw["pos"] = pd.to_numeric(chrom_pos[1], errors="coerce")
            distributions[group].append(overlap(draw.dropna(subset=["pos"]), intervals).mean())
    return {group: np.asarray(values) for group, values in distributions.items()}


def plot_group(group: str, observed: float, observed_n: int, observed_k: int, background: np.ndarray, outbase: Path) -> dict[str, float]:
    bg_mean = float(background.mean())
    bg_low, bg_high = np.quantile(background, [0.025, 0.975])
    obs_low, obs_high = wilson_interval(observed_k, observed_n)
    empirical_p = (1 + int(np.sum(background >= observed))) / (len(background) + 1)
    odds_ratio = (observed / (1 - observed)) / (bg_mean / (1 - bg_mean)) if 0 < observed < 1 and 0 < bg_mean < 1 else float("nan")
    fig, ax = plt.subplots(figsize=(3.15, 3.25))
    vals = np.array([observed, bg_mean]) * 100
    err = np.array([[observed - obs_low, bg_mean - bg_low], [obs_high - observed, bg_high - bg_mean]]) * 100
    bars = ax.bar([0, 1], vals, width=0.56, color=["#C75B5B", "#B8BEC5"], edgecolor="none", zorder=2)
    ax.errorbar([0, 1], vals, yerr=err, fmt="none", color="#222222", capsize=2.5, linewidth=0.8, zorder=3)
    ax.text(0, vals[0] + max(vals) * 0.04, f"{observed_k}/{observed_n}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks([0, 1], [group, "Matched cis\nbackground"])
    ax.set_ylabel("E5 enhancer overlap (%)")
    ax.set_ylim(0, max((vals + err[1]).max() * 1.28, 8))
    ax.tick_params(direction="out", length=3, width=0.7)
    ax.text(0, 1.085, f"{group} ieQTL", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=8.8, fontweight="bold")
    ax.text(0, 1.025, f"Observed {observed_k}/{observed_n}; OR = {odds_ratio:.2f}; P = {empirical_p:.3f}",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.8, color="#4A4A4A")
    fig.subplots_adjust(top=0.78, bottom=0.19, left=0.20, right=0.98)
    fig.savefig(outbase.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(outbase.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return {"group": group, "observed_overlap": observed_k, "observed_total": observed_n, "observed_rate": observed,
            "background_mean_rate": bg_mean, "background_low": bg_low, "background_high": bg_high,
            "odds_ratio": odds_ratio, "empirical_enrichment_p": empirical_p}


def main() -> None:
    args = parse_args()
    out = Path(args.out_dir)
    if out.exists():
        raise FileExistsError(f"Output directory already exists; choose a new path: {out}")
    leads_path = Path(args.leads)
    initial = pd.read_csv(leads_path, sep="\t")
    required = {"group", "tissue", "cell_type", "phenotype2", "variant_id"}
    missing = required.difference(initial.columns)
    if missing:
        raise ValueError(f"Lead table is missing columns: {sorted(missing)}")
    if initial.empty:
        raise ValueError("Lead table has no rows")
    initial = initial[initial.group.isin(["PPH3-high", "PPH4-high"])].copy()
    if initial.empty:
        raise ValueError("Lead table contains no PPH3-high or PPH4-high rows")
    out.mkdir(parents=True)
    annotated, choices = attach_nominal_and_background(initial, Path(args.ieqtl_root))
    chrom_pos = annotated.variant_id.astype(str).str.split("_", n=1, expand=True)
    annotated["chrom"] = "chr" + chrom_pos[0].str.replace("chr", "", regex=False)
    annotated["pos"] = pd.to_numeric(chrom_pos[1], errors="coerce")
    intervals = load_intervals(Path(args.e5_bed))
    annotated["overlap_tissue_specific_E5"] = overlap(annotated.dropna(subset=["pos"]), intervals)
    annotated.to_csv(out / "pph_high_ieqtl_leads_with_matched_background.tsv", sep="\t", index=False)
    distributions = build_background_distribution(annotated, choices, intervals, args.permutations, args.seed)
    rows = []
    for group in ["PPH3-high", "PPH4-high"]:
        data = annotated[annotated.group == group]
        observed_k = int(data.overlap_tissue_specific_E5.sum())
        rows.append(plot_group(group, observed_k / len(data), len(data), observed_k, distributions[group], out / f"{group}_ieqtl_vs_matched_background_E5"))
        pd.DataFrame({"background_overlap_rate": distributions[group]}).to_csv(out / f"{group}_matched_background_permutations.tsv", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(out / "pph_high_ieqtl_matched_background_E5_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
