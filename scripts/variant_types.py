"""Classify VCF variants into SNPs vs. INDELs and export results/variant_types.csv."""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

from read_vcf import ROOT, iter_vcf_records


def export_variant_types(vcf_path: Path, output_csv: Path) -> pd.DataFrame:
    """Classify each VCF record as SNP or INDEL and write summary CSV."""
    rows = [
        {
            "CHROM": r.chrom,
            "POS": r.pos,
            "ID": r.var_id,
            "REF": r.ref,
            "ALT": r.alt,
            "TYPE": r.variant_type,
        }
        for r in iter_vcf_records(vcf_path)
    ]
    df = pd.DataFrame(rows, columns=["CHROM", "POS", "ID", "REF", "ALT", "TYPE"])
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
        default=ROOT / "results/variant_types.csv",
    )
    args = parser.parse_args(argv)
    df = export_variant_types(args.input, args.output)
    print(df["TYPE"].value_counts())
    print(f"Saved variant classifications -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
