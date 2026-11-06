"""Plot chromosome variant count and 1-Mb genomic density along Chromosome 22."""

from __future__ import annotations

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from read_vcf import ROOT, iter_vcf_records


def plot_chromosome_distribution(vcf_path: Path, output_png: Path) -> None:
    """Render two-panel Chromosome 22 summary and 1-Mb genomic density profile."""
    records = list(iter_vcf_records(vcf_path))
    positions_mb = np.array([r.pos / 1e6 for r in records])
    known_count = sum(1 for r in records if r.is_known_dbsnp)
    novel_count = len(records) - known_count

    sns.set_theme(style="whitegrid", context="talk")
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(15, 6), dpi=300, gridspec_kw={"width_ratios": [1, 1.65]}
    )

    # Panel A: Chromosome 22 Variant Count Breakdown (Known dbSNP vs Novel)
    categories = ["Total Chr22\nVariants", "Known dbSNP\n(rsID)", "Novel Loci\n(Unannotated)"]
    values = [len(records), known_count, novel_count]
    colors = ["#2b5c8f", "#2ca02c", "#e67e22"]

    bars = ax1.bar(
        categories,
        values,
        color=colors,
        edgecolor="black",
        linewidth=0.9,
        width=0.55,
        alpha=0.88,
    )
    for bar in bars:
        h = int(bar.get_height())
        pct = 100.0 * h / len(records) if records else 0.0
        ax1.annotate(
            f"{h:,}\n({pct:.1f}%)",
            xy=(bar.get_x() + bar.get_width() / 2, h),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )
    ax1.set_title("A. Chromosome 22 Cohort Variant Summary", fontsize=13, fontweight="bold", pad=12)
    ax1.set_ylabel("Variant Count", fontsize=11, fontweight="bold")
    ax1.set_ylim(0, len(records) * 1.24)

    # Panel B: 1-Mb Window Variant Density Profile along GRCh37/GRCh38 Chromosome 22
    bins_mb = np.arange(0, 53, 1)
    counts_per_mb, edges = np.histogram(positions_mb, bins=bins_mb)
    centers = 0.5 * (edges[:-1] + edges[1:])

    ax2.axvspan(
        0,
        16.05,
        color="#95a5a6",
        alpha=0.25,
        label="22p Acrocentric / Centromeric Gap (0–16.0 Mb)",
    )
    ax2.bar(
        centers,
        counts_per_mb,
        width=0.85,
        color="#2b5c8f",
        edgecolor="white",
        alpha=0.88,
        label="Variant Density (Variants / 1 Mb Window)",
    )
    euchromatic_mean = counts_per_mb[counts_per_mb > 0].mean()
    ax2.axhline(
        euchromatic_mean,
        color="#d9534f",
        linestyle="--",
        linewidth=2.0,
        label=f"Mean 22q Density ({euchromatic_mean:.0f} vars/Mb)",
    )
    ax2.set_title(
        "B. Genomic Variant Density Across Chromosome 22 (1 Mb Windows)",
        fontsize=13,
        fontweight="bold",
        pad=12,
    )
    ax2.set_xlabel("Genomic Coordinate on Chromosome 22 (Mb)", fontsize=11, fontweight="bold")
    ax2.set_ylabel("Variants per 1 Mb Window", fontsize=11, fontweight="bold")
    ax2.set_xlim(0, 52.5)
    ax2.legend(frameon=True, facecolor="white", fontsize=9.5, loc="upper left")

    sns.despine(top=True, right=True)
    fig.suptitle(
        "1000 Genomes Project (HG00098) — Chromosome 22 Variant Distribution",
        fontsize=15,
        fontweight="bold",
        y=1.01,
    )
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved chromosome distribution figure -> {output_png}")


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
        default=ROOT / "results/chromosome_distribution.png",
    )
    args = parser.parse_args(argv)
    plot_chromosome_distribution(args.input, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
