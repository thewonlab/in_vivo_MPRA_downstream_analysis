#!/usr/bin/env python3

import argparse
import json
import logging
import re
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd


LOGGER = logging.getLogger(__name__)


def load_count_file(filepath: Path) -> Counter:
    counts = Counter()
    with filepath.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            parts = line.strip().split(maxsplit=1)
            if len(parts) != 2:
                continue
            count_str, barcode = parts
            try:
                counts[barcode] += int(count_str)
            except ValueError:
                LOGGER.warning("Invalid count at %s:%d -> %r", filepath, line_no, line.strip())
    return counts


def find_sample_files(count_dir: Path, sample: str, strict_prefix: bool = True) -> List[Path]:
    candidates = sorted(p for p in count_dir.iterdir() if p.name.endswith("_count.txt"))
    if strict_prefix:
        pattern = re.compile(rf"^{re.escape(sample)}.*_count\.txt$")
        return [p for p in candidates if pattern.match(p.name)]
    return [p for p in candidates if sample in p.name]


def build_count_matrix(count_dir: Path, samples: Iterable[str], strict_prefix: bool = True) -> pd.DataFrame:
    count_df = pd.DataFrame()
    for sample in sorted(samples):
        sample_files = find_sample_files(count_dir, sample, strict_prefix=strict_prefix)
        if not sample_files:
            LOGGER.warning("No count files found for sample: %s", sample)

        sample_counts = Counter()
        for filepath in sample_files:
            sample_counts += load_count_file(filepath)

        count_df[sample] = pd.Series(sample_counts, dtype="int64")

    return count_df.fillna(0).astype(int)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge per-lane barcode count files into a sample-level count matrix"
    )
    parser.add_argument("--count-dir", required=True, help="Directory containing *_count.txt files")
    parser.add_argument("--fastq-dict", required=True, help="Path to fastq_dict.json")
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument(
        "--non-strict-sample-match",
        action="store_true",
        help="Use substring matching for sample names (default: strict prefix match)",
    )
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return parser.parse_args()


def read_fastq_dict(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"fastq_dict must be a JSON object: {path}")
    return data


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level), format="[%(levelname)s] %(message)s")

    count_dir = Path(args.count_dir)
    fastq_dict_path = Path(args.fastq_dict)
    output_path = Path(args.output)

    if not count_dir.exists():
        raise FileNotFoundError(f"Count directory not found: {count_dir}")
    if not fastq_dict_path.exists():
        raise FileNotFoundError(f"fastq_dict.json not found: {fastq_dict_path}")

    fastq_dict = read_fastq_dict(fastq_dict_path)
    count_df = build_count_matrix(
        count_dir=count_dir,
        samples=fastq_dict.keys(),
        strict_prefix=not args.non_strict_sample_match,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    count_df.to_csv(output_path, index_label="Barcode")
    LOGGER.info("Wrote merged count matrix: %s", output_path)


if __name__ == "__main__":
    main()
