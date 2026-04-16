#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pool replicate-level RNA/DNA counts for downstream MPRA analysis")
    parser.add_argument("--input", required=True, help="Input replicate-level count matrix CSV")
    parser.add_argument("--output", required=True, help="Output directory for pooled matrices")
    parser.add_argument("--prefix", required=True, help="Prefix for pooled output files")
    parser.add_argument("--pooling_scheme", required=True, help="Pooling scheme CSV")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    ct_mat = pd.read_csv(args.input, index_col=0)
    pooling_scheme = pd.read_csv(args.pooling_scheme, index_col=0)
    pooling_scheme = pooling_scheme.reset_index().rename(columns={"index": "Rep_name"})

    pooling_scheme_dict = (
        pooling_scheme[pooling_scheme["Pool"].str.contains("Pool")]
        .groupby("Pool")["Rep_name"]
        .apply(list)
        .to_dict()
    )

    pooled_rna = pd.DataFrame(index=ct_mat.index, columns=pooling_scheme_dict.keys())
    pooled_dna = pd.DataFrame(index=ct_mat.index, columns=pooling_scheme_dict.keys())

    for pool_name, rep_names in pooling_scheme_dict.items():
        rna_name = [f"{name}_CORT_R" for name in rep_names]
        dna_name = [f"{name}_CORT_D" for name in rep_names]
        pooled_rna[pool_name] = ct_mat[rna_name].sum(axis=1)
        pooled_dna[pool_name] = ct_mat[dna_name].sum(axis=1)

    pooled_meta = pd.DataFrame(index=pooled_rna.columns)
    pooled_meta["Virus"] = pooling_scheme.groupby("Pool")["Virus"].apply(lambda x: ",".join(x))
    pooled_meta["Sex"] = pooling_scheme.groupby("Pool")["Sex"].apply(lambda x: ",".join(x))
    pooled_meta["PND"] = pooling_scheme.groupby("Pool")["PND"].apply(lambda x: ",".join(x))
    pooled_meta["RIN"] = pooling_scheme.groupby("Pool")["RIN"].mean()
    pooled_meta["batch"] = pooling_scheme.groupby("Pool")["batch"].apply(lambda x: ",".join(x))
    pooled_meta["sample_list"] = pooling_scheme.groupby("Pool")["Rep_name"].apply(lambda x: ",".join(x))

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    pooled_rna.to_csv(output_dir / f"{args.prefix}_RNA.csv")
    pooled_dna.to_csv(output_dir / f"{args.prefix}_DNA.csv")
    pooled_meta.to_csv(output_dir / f"{args.prefix}_meta.csv")


if __name__ == "__main__":
    main()
