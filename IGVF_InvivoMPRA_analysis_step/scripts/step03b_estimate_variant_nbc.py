#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Estimate per-variant barcode counts from the v2 barcode mapping tables used to build count_matrix_variant_merged_v2_cortex.csv"
    )
    parser.add_argument("--barcode-mapping-vb", required=True, help="Path to barcode_mapping_vb_v2.csv")
    parser.add_argument("--barcode-mapping-gvvc2", required=True, help="Path to barcode_mapping_gvvc2_v2.csv")
    parser.add_argument(
        "--variant-matrix",
        help="Optional path to count_matrix_variant_merged_v2_cortex.csv; if provided, output is restricted to variants present there",
    )
    parser.add_argument("--output", required=True, help="Output CSV path with columns: variant,n_bc")
    return parser.parse_args()


def load_variant_set(path: Path) -> pd.Index:
    df = pd.read_csv(path, usecols=["Variant"])
    return pd.Index(df["Variant"].astype(str))


def main() -> None:
    args = parse_args()

    vb = pd.read_csv(args.barcode_mapping_vb, usecols=["variant"])
    gvvc2 = pd.read_csv(args.barcode_mapping_gvvc2, usecols=["variant"])

    counts = (
        pd.concat([vb["variant"], gvvc2["variant"]], axis=0)
        .dropna()
        .astype(str)
        .value_counts()
        .rename_axis("variant")
        .reset_index(name="n_bc")
    )

    if args.variant_matrix:
        keep = load_variant_set(Path(args.variant_matrix))
        counts = keep.to_frame(index=False, name="variant").merge(counts, on="variant", how="left")
        counts["n_bc"] = counts["n_bc"].fillna(0).astype(int)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_path, index=False)


if __name__ == "__main__":
    main()
