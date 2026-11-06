"""VCF v4.x parser and variant record extractor for germline/cohort QC analysis."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

ROOT = Path(__file__).resolve().parents[1]
TRANSITION_PAIRS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


@dataclass(frozen=True)
class VcfRecord:
    chrom: str
    pos: int
    var_id: str
    ref: str
    alt: str
    qual: str
    filter_status: str
    genotype: str

    @property
    def is_known_dbsnp(self) -> bool:
        return self.var_id.startswith("rs")

    @property
    def variant_type(self) -> str:
        if len(self.ref) == 1 and len(self.alt) == 1:
            return "SNP"
        return "INDEL"

    @property
    def substitution_class(self) -> str:
        if self.variant_type != "SNP":
            return "INDEL"
        return "Transition" if (self.ref, self.alt) in TRANSITION_PAIRS else "Transversion"

    @property
    def zygosity(self) -> str:
        gt = self.genotype
        if gt in {"0|1", "1|0", "0/1", "1/0", "2|1", "1|2"}:
            return "Heterozygous"
        if gt in {"1|1", "1/1", "2|2", "2/2"}:
            return "Homozygous_Alt"
        if gt in {"0|0", "0/0"}:
            return "Homozygous_Ref"
        return "Other"


def iter_vcf_records(vcf_path: Path) -> Iterator[VcfRecord]:
    """Stream non-header variant records from a VCF file."""
    with vcf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 8:
                continue
            gt = parts[9].split(":")[0] if len(parts) > 9 else "."
            yield VcfRecord(
                chrom=parts[0],
                pos=int(parts[1]),
                var_id=parts[2],
                ref=parts[3],
                alt=parts[4],
                qual=parts[5],
                filter_status=parts[6],
                genotype=gt,
            )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "data/raw/sample.vcf",
    )
    args = parser.parse_args(argv)
    records = list(iter_vcf_records(args.input))
    known = sum(1 for r in records if r.is_known_dbsnp)
    print(f"Total variants parsed : {len(records):,}")
    if records:
        first = records[0]
        print(
            f"First variant         : chr{first.chrom}:{first.pos} {first.var_id} "
            f"{first.ref}>{first.alt} ({first.genotype})"
        )
        print(f"Known dbSNP (rsID)    : {known:,} ({100.0 * known / len(records):.2f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
