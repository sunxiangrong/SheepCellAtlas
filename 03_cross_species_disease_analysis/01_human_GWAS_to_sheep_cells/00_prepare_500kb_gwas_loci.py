#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def find_mlma(trait_dir: Path) -> Path | None:
    files = sorted(trait_dir.glob("*_ALL_MERGED_LOCI.mlma"))
    return files[0] if files else None


def crop_trait(in_file: Path, out_file: Path, window: int, cluster_gap: int) -> dict:
    df = pd.read_csv(in_file, sep=r"\s+|\t|,", engine="python")
    required = {"Chr", "SNP", "bp", "p"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{in_file} missing columns: {sorted(missing)}")

    df = df.copy()
    df["bp"] = pd.to_numeric(df["bp"], errors="coerce")
    df["p"] = pd.to_numeric(df["p"], errors="coerce")
    df = df.dropna(subset=["Chr", "SNP", "bp", "p"])
    df["bp"] = df["bp"].astype(int)

    kept_parts = []
    locus_count = 0
    for chrom, sub in df.groupby(df["Chr"].astype(str), sort=False):
        sub = sub.sort_values("bp").copy()
        if sub.empty:
            continue
        cluster_id = (sub["bp"].diff().fillna(0) > cluster_gap).cumsum()
        for _, block in sub.groupby(cluster_id, sort=False):
            if block.empty:
                continue
            locus_count += 1
            lead = block.loc[block["p"].idxmin()]
            lo = int(lead["bp"]) - window
            hi = int(lead["bp"]) + window
            kept_parts.append(block[(block["bp"] >= lo) & (block["bp"] <= hi)].copy())

    if kept_parts:
        out = pd.concat(kept_parts, ignore_index=True).drop_duplicates(subset=["Chr", "SNP", "bp"])
        out = out.sort_values(["Chr", "bp"], key=lambda s: s.astype(str) if s.name == "Chr" else s)
    else:
        out = df.iloc[0:0].copy()

    out_file.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_file, sep="\t", index=False)
    return {
        "trait": in_file.parent.name,
        "input_file": str(in_file),
        "output_file": str(out_file),
        "input_rows": int(len(df)),
        "output_rows": int(len(out)),
        "locus_clusters": int(locus_count),
        "window_bp_each_side": int(window),
        "cluster_gap_bp": int(cluster_gap),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", required=True)
    parser.add_argument("--out-root", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--window", type=int, default=500_000)
    parser.add_argument("--cluster-gap", type=int, default=2_500_000)
    args = parser.parse_args()

    input_root = Path(args.input_root)
    out_root = Path(args.out_root)
    rows = []
    for trait_dir in sorted([p for p in input_root.iterdir() if p.is_dir()]):
        in_file = find_mlma(trait_dir)
        if in_file is None:
            rows.append({"trait": trait_dir.name, "status": "missing_mlma"})
            continue
        out_file = out_root / trait_dir.name / f"{trait_dir.name}_ALL_MERGED_LOCI.mlma"
        row = crop_trait(in_file, out_file, args.window, args.cluster_gap)
        row["status"] = "ok" if row["output_rows"] > 0 else "empty_after_crop"
        rows.append(row)

    manifest = pd.DataFrame(rows)
    Path(args.manifest).parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.manifest, sep="\t", index=False)
    print(f"traits={len(manifest)} ok={(manifest['status'] == 'ok').sum()} out={out_root}")


if __name__ == "__main__":
    main()
