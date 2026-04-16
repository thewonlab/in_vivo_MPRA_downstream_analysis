#!/usr/bin/env python3

import argparse
import gzip
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

GRCH38_ACCESSIONS = {
    "chr1": "NC_000001.11",
    "chr2": "NC_000002.12",
    "chr3": "NC_000003.12",
    "chr4": "NC_000004.12",
    "chr5": "NC_000005.10",
    "chr6": "NC_000006.12",
    "chr7": "NC_000007.14",
    "chr8": "NC_000008.11",
    "chr9": "NC_000009.12",
    "chr10": "NC_000010.11",
    "chr11": "NC_000011.10",
    "chr12": "NC_000012.12",
    "chr13": "NC_000013.11",
    "chr14": "NC_000014.9",
    "chr15": "NC_000015.10",
    "chr16": "NC_000016.10",
    "chr17": "NC_000017.11",
    "chr18": "NC_000018.10",
    "chr19": "NC_000019.10",
    "chr20": "NC_000020.11",
    "chr21": "NC_000021.9",
    "chr22": "NC_000022.11",
    "chrX": "NC_000023.11",
    "chrY": "NC_000024.10",
    "chrM": "NC_012920.1",
    "chrMT": "NC_012920.1",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Format pooled MPRA outputs into IGVF submission tables"
    )
    parser.add_argument("--config-json", required=True, help="Path to step6 config JSON")
    return parser.parse_args()


def read_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def load_sequence_design(path: Path | None) -> pd.DataFrame | None:
    if path is None or not path.exists():
        return None

    df = pd.read_csv(path, sep="\t", dtype=str)
    for col in ["start", "end"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    if "variant_pos" in df.columns:
        df["variant_pos"] = (
            df["variant_pos"]
            .fillna("")
            .astype(str)
            .str.replace("[", "", regex=False)
            .str.replace("]", "", regex=False)
            .replace("", pd.NA)
        )
        df["variant_pos"] = pd.to_numeric(df["variant_pos"], errors="coerce").astype("Int64")
    return df


def load_variant_nbc(path: Path | None) -> pd.Series | None:
    if path is None or not path.exists():
        return None
    df = pd.read_csv(path)
    if not {"variant", "n_bc"}.issubset(df.columns):
        return None
    series = pd.to_numeric(df.set_index("variant")["n_bc"], errors="coerce")
    series = series.where(series > 0, np.nan)
    return series


def apply_name_suffix(filename: str, suffix: str) -> str:
    if not suffix:
        return filename
    if filename.endswith(".tsv.gz"):
        return filename[:-7] + suffix + ".tsv.gz"
    if filename.endswith(".bed.gz"):
        return filename[:-7] + suffix + ".bed.gz"
    return filename + suffix


def normalize_cpm(df: pd.DataFrame) -> pd.DataFrame:
    column_sums = df.sum(axis=0)
    return df.divide(column_sums, axis=1) * 1_000_000


def neg_log10(series: pd.Series) -> pd.Series:
    series = pd.to_numeric(series, errors="coerce")
    with np.errstate(divide="ignore"):
        out = -np.log10(series)
    return pd.Series(out, index=series.index)


def parse_variant_id(variant_id: str) -> dict:
    variant_id = str(variant_id)
    parts = variant_id.split("_", 2)
    chrom = parts[0] if len(parts) > 0 else ""
    pos = np.nan
    if len(parts) > 1:
        try:
            pos = int(parts[1])
        except ValueError:
            pos = np.nan

    ref_allele = ""
    alt_allele = ""
    if ":" in variant_id:
        tail = variant_id.rsplit("_", 1)[-1]
        allele_parts = tail.split(":")
        if len(allele_parts) >= 2:
            ref_allele = allele_parts[-2]
            alt_allele = allele_parts[-1]

    return {
        "chrom": chrom,
        "pos": pos,
        "refAllele": ref_allele,
        "altAllele": alt_allele,
    }


def build_grch38_spdi(chrom: str, start: int, variant_pos: int, ref_allele: str, alt_allele: str) -> str | None:
    accession = GRCH38_ACCESSIONS.get(str(chrom))
    if accession is None:
        return None
    if pd.isna(start) or pd.isna(variant_pos):
        return None
    if not ref_allele or not alt_allele:
        return None
    pos0 = int(start) + int(variant_pos)
    return f"{accession}:{pos0}:{ref_allele}:{alt_allele}"


def build_variant_metadata(sequence_design: pd.DataFrame | None, reporter_variant: pd.DataFrame) -> pd.DataFrame:
    variant_col = "variant_id" if "variant_id" in reporter_variant.columns else "variant_name"
    base = reporter_variant[[variant_col, "refAllele", "altAllele"]].copy()
    base = base.rename(columns={variant_col: "variant_id"})
    base["spdi"] = pd.NA
    base["chrom"] = pd.NA
    base["chromStart"] = pd.NA
    base["chromEnd"] = pd.NA
    base["strand"] = pd.NA

    if sequence_design is None:
        return base

    seq = sequence_design.copy()
    seq["variant_id"] = seq["name"].astype(str).str.rsplit("_", n=1).str[0]
    seq["variant_pos_clean"] = seq["variant_pos"]

    alt_rows = seq[seq["allele"] == "alt"].copy()
    alt_rows["spdi_built"] = alt_rows.apply(
        lambda row: build_grch38_spdi(
            chrom=row.get("chr"),
            start=row.get("start"),
            variant_pos=row.get("variant_pos_clean"),
            ref_allele=row.get("SPDI", "").split(":")[-2] if pd.notna(row.get("SPDI")) and len(str(row.get("SPDI")).split(":")) >= 4 else "",
            alt_allele=row.get("SPDI", "").split(":")[-1] if pd.notna(row.get("SPDI")) and len(str(row.get("SPDI")).split(":")) >= 4 else "",
        ),
        axis=1,
    )
    alt_rows["spdi"] = alt_rows["SPDI"].where(alt_rows["SPDI"].notna(), alt_rows["spdi_built"])
    alt_rows = alt_rows.rename(columns={"chr": "chrom", "start": "chromStart", "end": "chromEnd"})

    meta = (
        alt_rows[["variant_id", "spdi", "chrom", "chromStart", "chromEnd", "strand"]]
        .drop_duplicates(subset=["variant_id"])
        .copy()
    )

    out = base.merge(meta, on="variant_id", how="left", suffixes=("", "_meta"))
    for col in ["spdi", "chrom", "chromStart", "chromEnd", "strand"]:
        meta_col = f"{col}_meta"
        if meta_col in out.columns:
            out[col] = out[meta_col].combine_first(out[col])
            out = out.drop(columns=[meta_col])

    missing = out["spdi"].isna()
    if missing.any():
        parsed = out.loc[missing, "variant_id"].map(parse_variant_id).apply(pd.Series)
        out.loc[missing, "chrom"] = parsed["chrom"].values
        out.loc[missing, "chromStart"] = (parsed["pos"] - 75).astype("Int64").values
        out.loc[missing, "chromEnd"] = (parsed["pos"] + 75).astype("Int64").values
        out.loc[missing, "strand"] = "+"
        out.loc[missing, "spdi"] = out.loc[missing].apply(
            lambda row: build_grch38_spdi(
                chrom=row["chrom"],
                start=row["chromStart"],
                variant_pos=75,
                ref_allele=row["refAllele"],
                alt_allele=row["altAllele"],
            ),
            axis=1,
        )

    return out


def build_reporter_experiment(
    rna: pd.DataFrame,
    dna: pd.DataFrame,
    variant_nbc: pd.Series | None = None,
) -> pd.DataFrame:
    dna_cpm = normalize_cpm(dna)
    rna_cpm = normalize_cpm(rna)

    experiment_rows = []
    for replicate in dna.columns:
        frame = pd.DataFrame(
            {
                "replicate": replicate,
                "oligo_name": dna.index,
                "dna_counts": dna[replicate].values,
                "rna_counts": rna[replicate].values,
                "dna_normalized": dna_cpm[replicate].values,
                "rna_normalized": rna_cpm[replicate].values,
            }
        )
        frame["log2FoldChange"] = np.log2(frame["rna_normalized"] / frame["dna_normalized"])
        if variant_nbc is not None:
            frame["n_bc"] = frame["oligo_name"].map(variant_nbc)
        else:
            frame["n_bc"] = np.nan
        experiment_rows.append(frame)

    return pd.concat(experiment_rows, ignore_index=True)


def build_reporter_element(
    mpra_element: pd.DataFrame, rna_cpm: pd.DataFrame, dna_cpm: pd.DataFrame
) -> pd.DataFrame:
    out = mpra_element.copy()
    out["oligo_name"] = out["target"]
    out["log2FoldChange"] = out["logFC_active"]
    out["inputCount"] = out["target"].map(dna_cpm.mean(axis=1))
    out["outputCount"] = out["target"].map(rna_cpm.mean(axis=1))
    out["minusLog10PValue"] = neg_log10(out["pVal_active"])
    out["minusLog10QValue"] = neg_log10(out["FDR_active"])
    return out[
        [
            "oligo_name",
            "log2FoldChange",
            "inputCount",
            "outputCount",
            "minusLog10PValue",
            "minusLog10QValue",
        ]
    ]


def build_reporter_variant(
    mpra_variant: pd.DataFrame,
    rna_cpm: pd.DataFrame,
    dna_cpm: pd.DataFrame,
    sequence_design: pd.DataFrame | None = None,
) -> pd.DataFrame:
    variant_id_col = "Unnamed: 0" if "Unnamed: 0" in mpra_variant.columns else mpra_variant.columns[0]
    out = mpra_variant.copy()
    out["variant_name"] = out[variant_id_col]

    parsed = out["variant_name"].map(parse_variant_id).apply(pd.Series)
    out["refAllele"] = parsed["refAllele"]
    out["altAllele"] = parsed["altAllele"]

    out["ref_oligo"] = out["variant_name"] + "_" + out["refAllele"]
    out["alt_oligo"] = out["variant_name"] + "_" + out["altAllele"]

    dna_mean = dna_cpm.mean(axis=1)
    rna_mean = rna_cpm.mean(axis=1)

    out["log2FoldChange"] = out["logFC_allelic"]
    out["inputCountRef"] = out["ref_oligo"].map(dna_mean)
    out["outputCountRef"] = out["ref_oligo"].map(rna_mean)
    out["inputCountAlt"] = out["alt_oligo"].map(dna_mean)
    out["outputCountAlt"] = out["alt_oligo"].map(rna_mean)
    out["minusLog10PValue"] = neg_log10(out["pVal_allelic"])
    out["minusLog10QValue"] = neg_log10(out["FDR_allelic"])
    if "postProbEffect" in out.columns:
        out["postProbEffect"] = pd.to_numeric(out["postProbEffect"], errors="coerce")
    elif "B_allelic" in out.columns:
        out["postProbEffect"] = 1 / (1 + np.exp(-pd.to_numeric(out["B_allelic"], errors="coerce")))
    else:
        out["postProbEffect"] = np.nan

    if {"CI_lower_95", "CI_upper_95"}.issubset(out.columns):
        out["CI_lower_95"] = pd.to_numeric(out["CI_lower_95"], errors="coerce")
        out["CI_upper_95"] = pd.to_numeric(out["CI_upper_95"], errors="coerce")
    elif {"logFC_allelic", "t_allelic"}.issubset(out.columns):
        logfc = pd.to_numeric(out["logFC_allelic"], errors="coerce")
        t_stat = pd.to_numeric(out["t_allelic"], errors="coerce")
        se = (logfc.abs() / t_stat.abs()).replace([np.inf, -np.inf], np.nan)
        out["CI_lower_95"] = logfc - 1.96 * se
        out["CI_upper_95"] = logfc + 1.96 * se
    else:
        out["CI_lower_95"] = np.nan
        out["CI_upper_95"] = np.nan
    out["variantPos"] = 75

    variant_meta = build_variant_metadata(sequence_design, out[["variant_name", "refAllele", "altAllele"]].rename(columns={"variant_name": "variant_id"}))
    variant_meta = variant_meta.rename(columns={"variant_id": "variant_name"})
    out = out.merge(variant_meta[["variant_name", "spdi"]], on="variant_name", how="left")
    out["variant_id"] = out["spdi"].combine_first(out["variant_name"])

    return out[
        [
            "variant_id",
            "log2FoldChange",
            "inputCountRef",
            "outputCountRef",
            "inputCountAlt",
            "outputCountAlt",
            "minusLog10PValue",
            "minusLog10QValue",
            "postProbEffect",
            "CI_lower_95",
            "CI_upper_95",
            "variantPos",
            "refAllele",
            "altAllele",
        ]
    ]


def build_genomic_element_effect(
    reporter_element: pd.DataFrame, sequence_design: pd.DataFrame | None = None
) -> pd.DataFrame:
    out = reporter_element.copy()
    if sequence_design is not None:
        seq = sequence_design.rename(
            columns={"name": "oligo_name", "chr": "chrom", "start": "chromStart", "end": "chromEnd"}
        )
        out = out.merge(
            seq[["oligo_name", "chrom", "chromStart", "chromEnd", "strand"]],
            on="oligo_name",
            how="left",
        )
    else:
        parsed = out["oligo_name"].str.rsplit("_", n=1).str[0].map(parse_variant_id).apply(pd.Series)
        out["chrom"] = parsed["chrom"]
        out["chromStart"] = (parsed["pos"] - 75).astype("Int64")
        out["chromEnd"] = (parsed["pos"] + 75).astype("Int64")
        out["strand"] = "."

    out["name"] = out["oligo_name"]
    out["score"] = 0
    return out[
        [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "log2FoldChange",
            "inputCount",
            "outputCount",
            "minusLog10PValue",
            "minusLog10QValue",
        ]
    ]


def build_genomic_variant_effect(
    reporter_variant: pd.DataFrame, sequence_design: pd.DataFrame | None = None
) -> pd.DataFrame:
    out = reporter_variant.copy()
    if sequence_design is not None:
        seq = sequence_design.copy()
        seq = seq[seq["allele"] == "alt"].copy()
        seq["spdi"] = seq["SPDI"]
        seq = seq.rename(columns={"chr": "chrom", "start": "chromStart", "end": "chromEnd"})
        seq = seq[["spdi", "chrom", "chromStart", "chromEnd", "strand"]].drop_duplicates(subset=["spdi"])
        out = out.merge(
            seq,
            left_on="variant_id",
            right_on="spdi",
            how="left",
        )
        out = out.drop(columns=["spdi"])
    else:
        parsed = out["variant_id"].map(parse_variant_id).apply(pd.Series)
        out["chrom"] = parsed["chrom"]
        out["chromStart"] = (parsed["pos"] - 75).astype("Int64")
        out["chromEnd"] = (parsed["pos"] + 75).astype("Int64")
        out["strand"] = "."

    out["name"] = out["variant_id"]
    out["score"] = 0
    return out[
        [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "log2FoldChange",
            "inputCountRef",
            "outputCountRef",
            "inputCountAlt",
            "outputCountAlt",
            "minusLog10PValue",
            "minusLog10QValue",
            "postProbEffect",
            "CI_lower_95",
            "CI_upper_95",
            "variantPos",
            "refAllele",
            "altAllele",
        ]
    ]


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip")


def write_bed_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        df.to_csv(handle, sep="\t", index=False, header=False)


def maybe_liftover_bed_gz(
    input_bed_gz: Path,
    output_bed_gz: Path,
    liftover_executable: str | None,
    chain_file: str | None,
) -> bool:
    if not liftover_executable or not chain_file:
        return False

    liftover_path = shutil.which(liftover_executable) or liftover_executable
    if not Path(liftover_path).exists():
        return False

    chain_path = Path(chain_file)
    if not chain_path.exists():
        return False

    output_bed_gz.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        input_bed = tmpdir / "input.bed"
        lifted_bed = tmpdir / "lifted.bed"
        unmapped_bed = tmpdir / "unmapped.bed"

        with gzip.open(input_bed_gz, "rt", encoding="utf-8") as src, input_bed.open("w", encoding="utf-8") as dst:
            dst.write(src.read())

        subprocess.run(
            [str(liftover_path), str(input_bed), str(chain_path), str(lifted_bed), str(unmapped_bed)],
            check=True,
        )

        with lifted_bed.open("r", encoding="utf-8") as src, gzip.open(output_bed_gz, "wt", encoding="utf-8") as dst:
            dst.write(src.read())

    return True


def main() -> None:
    args = parse_args()
    config = read_config(Path(args.config_json))

    pooled_rna = read_matrix(Path(config["pooled_rna"]))
    pooled_dna = read_matrix(Path(config["pooled_dna"]))
    mpra_element = pd.read_csv(config["mpra_element"])
    mpra_variant = pd.read_csv(config["mpra_variant"])
    sequence_design = load_sequence_design(
        Path(config["sequence_design_tsv"]) if config.get("sequence_design_tsv") else None
    )
    variant_nbc = load_variant_nbc(
        Path(config["variant_nbc_csv"]) if config.get("variant_nbc_csv") else None
    )
    output_dir = Path(config["output_dir"])
    output_name_suffix = str(config.get("output_name_suffix", ""))

    rna_cpm = normalize_cpm(pooled_rna)
    dna_cpm = normalize_cpm(pooled_dna)

    reporter_experiment = build_reporter_experiment(pooled_rna, pooled_dna, variant_nbc=variant_nbc)
    reporter_element = build_reporter_element(mpra_element, rna_cpm, dna_cpm)
    reporter_variant = build_reporter_variant(mpra_variant, rna_cpm, dna_cpm, sequence_design=sequence_design)
    genomic_element = build_genomic_element_effect(reporter_element, sequence_design=sequence_design)
    genomic_variant = build_genomic_variant_effect(reporter_variant, sequence_design=sequence_design)

    write_tsv_gz(reporter_experiment, output_dir / apply_name_suffix("Reporter_Experiment.tsv.gz", output_name_suffix))
    write_tsv_gz(reporter_element, output_dir / apply_name_suffix("Reporter_Element.tsv.gz", output_name_suffix))
    write_tsv_gz(reporter_variant, output_dir / apply_name_suffix("Reporter_Variant.tsv.gz", output_name_suffix))
    write_bed_gz(genomic_element, output_dir / apply_name_suffix("Reporter_Genomic_Element_Effect.bed.gz", output_name_suffix))
    write_bed_gz(genomic_variant, output_dir / apply_name_suffix("Reporter_Genomic_Variant_Effect.bed.gz", output_name_suffix))

    if sequence_design is not None:
        write_bed_gz(genomic_element, output_dir / apply_name_suffix("Reporter_Genomic_Element_Effect.hg38.bed.gz", output_name_suffix))
        write_bed_gz(genomic_variant, output_dir / apply_name_suffix("Reporter_Genomic_Variant_Effect.hg38.bed.gz", output_name_suffix))

    maybe_liftover_bed_gz(
        input_bed_gz=output_dir / apply_name_suffix("Reporter_Genomic_Element_Effect.bed.gz", output_name_suffix),
        output_bed_gz=output_dir / apply_name_suffix("Reporter_Genomic_Element_Effect.hg38.bed.gz", output_name_suffix),
        liftover_executable=config.get("liftover_executable"),
        chain_file=config.get("hg19_to_hg38_chain"),
    )
    maybe_liftover_bed_gz(
        input_bed_gz=output_dir / apply_name_suffix("Reporter_Genomic_Variant_Effect.bed.gz", output_name_suffix),
        output_bed_gz=output_dir / apply_name_suffix("Reporter_Genomic_Variant_Effect.hg38.bed.gz", output_name_suffix),
        liftover_executable=config.get("liftover_executable"),
        chain_file=config.get("hg19_to_hg38_chain"),
    )


if __name__ == "__main__":
    main()
