"""Export structured variant table (vcf_summary.csv) and print cohort QC statistics."""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

from read_vcf import ROOT, iter_vcf_records


def generate_vcf_summary(vcf_path: Path, output_csv: Path) -> pd.DataFrame:
    """Parse VCF records into a canonical summary DataFrame and save to CSV."""
    rows = [
        {
            "CHROM": r.chrom,
            "POS": r.pos,
            "ID": r.var_id,
            "REF": r.ref,
            "ALT": r.alt,
        }
        for r in iter_vcf_records(vcf_path)
    ]
    df = pd.DataFrame(rows, columns=["CHROM", "POS", "ID", "REF", "ALT"])
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    return df


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "data/raw/sample.vcf",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results/vcf_summary.csv",
    )
    args = parser.parse_args(argv)

    df = generate_vcf_summary(args.input, args.output)
    print(df.head())
    print()
    print(df["CHROM"].value_counts())
    print(f"\nExported {len(df):,} variant records -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
