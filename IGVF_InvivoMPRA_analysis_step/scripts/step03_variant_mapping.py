#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def rev_comp(seq: str) -> str:
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(complement.get(base, base) for base in reversed(str(seq)))


def filtering_bc_var(df: pd.DataFrame, variant_col: str, min_obs_per_var: int = 5, min_nonmissing_frac: float = 0.5) -> pd.DataFrame:
    dna_col = [col for col in df.columns if col.endswith("_D")]
    ncol = (df.shape[1] - 1) / 2

    df_temp = df.drop(columns=[variant_col])
    bc_stat = df_temp[dna_col].notna().sum(axis=1)
    df_temp = df_temp.loc[bc_stat > ncol * min_nonmissing_frac].copy()
    df_temp[variant_col] = df.loc[df_temp.index, variant_col]

    var_counts = df_temp[variant_col].value_counts()
    keep_vars = var_counts[var_counts > min_obs_per_var].index
    return df_temp[df_temp[variant_col].isin(keep_vars)].copy()


def load_mapping(path: Path) -> dict:
    mapping = pd.read_csv(path, index_col=0)
    mapping.index = mapping.index.map(lambda x: rev_comp(str(x)))
    return mapping.to_dict()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Annotate barcode count matrix with variant IDs and build variant-level matrices")
    parser.add_argument("--count-matrix-merged", required=True, help="Merged barcode count matrix CSV")
    parser.add_argument("--barcode-mapping-vb", required=True, help="Barcode mapping CSV for VB")
    parser.add_argument("--barcode-mapping-gvvc1", required=True, help="Barcode mapping CSV for GVVC1")
    parser.add_argument("--barcode-mapping-gvvc2", required=True, help="Barcode mapping CSV for GVVC2")
    parser.add_argument("--output-with-variant", required=True, help="Annotated barcode-level output CSV")
    parser.add_argument("--output-variant-full", required=True, help="Variant-level full output CSV")
    parser.add_argument("--output-variant-filtered", required=True, help="Variant-level filtered output CSV")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    ct_merged = pd.read_csv(args.count_matrix_merged, index_col=0)

    vb_col = [f"VB{ind}_CORT" for ind in [1, 2, 3, 4, 5]]
    vb_col2 = [f"IHP{ind}_CORT" for ind in [3, 4, 7, 8, 9]]
    vb_col3 = [f"IHP{ind}_LIV" for ind in [6, 7, 8, 9]]
    gvvc1_col = [f"GVVC{ind}_CORT" for ind in [1, 3, 4, 5]]
    gvvc2_col = [f"GVVC{ind}_CORT" for ind in range(6, 23)]

    vb_bc_dict = load_mapping(Path(args.barcode_mapping_vb))
    gvvc1_bc_dict = load_mapping(Path(args.barcode_mapping_gvvc1))
    gvvc2_bc_dict = load_mapping(Path(args.barcode_mapping_gvvc2))

    ct_merged["variant_vb"] = ct_merged.index.map(vb_bc_dict["variant"])
    ct_merged["variant_gvvc1"] = ct_merged.index.map(gvvc1_bc_dict["variant"])
    ct_merged["variant_gvvc2"] = ct_merged.index.map(gvvc2_bc_dict["variant"])

    variant_index = list(
        set(ct_merged["variant_vb"].dropna().unique()).union(
            ct_merged["variant_gvvc1"].dropna().unique(),
            ct_merged["variant_gvvc2"].dropna().unique(),
        )
    )

    Path(args.output_with_variant).parent.mkdir(parents=True, exist_ok=True)
    ct_merged.to_csv(args.output_with_variant)

    target_col = [vb_col, gvvc1_col, gvvc2_col, vb_col2, vb_col3]
    target_virus = ["variant_vb", "variant_gvvc1", "variant_gvvc2", "variant_vb", "variant_vb"]

    merged_var_mat = pd.DataFrame(index=variant_index)
    for col_list, virus in zip(target_col, target_virus):
        sample_cols = [f"{c}_R" for c in col_list] + [f"{c}_D" for c in col_list] + [virus]
        existing_cols = [col for col in sample_cols if col in ct_merged.columns]
        if virus not in existing_cols:
            continue

        ct_temp = ct_merged[existing_cols].copy()
        ct_temp_var = filtering_bc_var(ct_temp, virus).groupby(virus).sum()
        merged_var_mat = merged_var_mat.join(ct_temp_var, how="outer")

    merged_var_mat_filtered = merged_var_mat.dropna()
    Path(args.output_variant_full).parent.mkdir(parents=True, exist_ok=True)
    merged_var_mat.to_csv(args.output_variant_full)
    merged_var_mat_filtered.to_csv(args.output_variant_filtered)


if __name__ == "__main__":
    main()
