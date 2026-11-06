"""Calculate Transition/Transversion (Ti/Tv) ratio and substitution spectrum from VCF."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from read_vcf import ROOT, TRANSITION_PAIRS, iter_vcf_records


def compute_titv_metrics(vcf_path: Path) -> dict[str, object]:
    """Compute Ti, Tv, Ti/Tv ratio, and 12-substitution counts across biallelic SNPs."""
    ti = 0
    tv = 0
    pair_counts: Counter[tuple[str, str]] = Counter()

    for rec in iter_vcf_records(vcf_path):
        if rec.variant_type == "SNP":
            pair = (rec.ref, rec.alt)
            pair_counts[pair] += 1
            if pair in TRANSITION_PAIRS:
                ti += 1
            else:
                tv += 1

    ratio = ti / tv if tv != 0 else 0.0
    return {
        "transitions": ti,
        "transversions": tv,
        "titv_ratio": ratio,
        "pair_counts": pair_counts,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "data/raw/sample.vcf",
    )
    args = parser.parse_args(argv)
    metrics = compute_titv_metrics(args.input)
    ti = metrics["transitions"]
    tv = metrics["transversions"]
    ratio = metrics["titv_ratio"]
    print(f"Transitions   : {ti:,} ({100.0 * ti / (ti + tv):.2f}%)")
    print(f"Transversions : {tv:,} ({100.0 * tv / (ti + tv):.2f}%)")
    print(f"Ti/Tv ratio   : {ratio:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
