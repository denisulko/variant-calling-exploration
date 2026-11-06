"""Plot variant type classification (SNP vs INDEL) and sample genotype zygosity."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

from read_vcf import ROOT, iter_vcf_records


def plot_variant_types(vcf_path: Path, output_png: Path) -> None:
    """Render two-panel figure of Variant Classes and Sample Zygosity (Het vs Hom-Alt)."""
    records = list(iter_vcf_records(vcf_path))
    total = len(records)

    class_counts = Counter()
    zygosity_counts = Counter()
    for r in records:
        if r.variant_type == "SNP":
            key = "Known dbSNP SNP" if r.is_known_dbsnp else "Novel SNP"
        else:
            key = "Multi-Allelic / INDEL"
        class_counts[key] += 1
        zygosity_counts[r.zygosity] += 1

    sns.set_theme(style="whitegrid", context="talk")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 6), dpi=300)

    # Panel A: Variant Type & dbSNP Annotation
    labels_a = ["Known dbSNP SNP", "Novel SNP", "Multi-Allelic / INDEL"]
    vals_a = [class_counts[k] for k in labels_a]
    colors_a = ["#2b5c8f", "#17becf", "#d9534f"]

    bars_a = ax1.bar(
        labels_a,
        vals_a,
        color=colors_a,
        edgecolor="black",
        linewidth=0.9,
        width=0.52,
        alpha=0.88,
    )
    for bar in bars_a:
        h = int(bar.get_height())
        pct = 100.0 * h / total if total else 0.0
        ax1.annotate(
            f"{h:,}\n({pct:.2f}%)",
            xy=(bar.get_x() + bar.get_width() / 2, h),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )
    ax1.set_title("A. Variant Classification & dbSNP Status", fontsize=13, fontweight="bold", pad=12)
    ax1.set_ylabel("Variant Count", fontsize=11, fontweight="bold")
    ax1.set_ylim(0, max(vals_a) * 1.24)

    # Panel B: Sample Genotype Zygosity (HG00098)
    het = zygosity_counts["Heterozygous"]
    hom_alt = zygosity_counts["Homozygous_Alt"]
    het_hom_ratio = het / hom_alt if hom_alt else 0.0

    labels_b = ["Heterozygous\n(0|1, 1|0)", "Homozygous Alt\n(1|1)"]
    vals_b = [het, hom_alt]
    colors_b = ["#2ca02c", "#8e44ad"]

    bars_b = ax2.bar(
        labels_b,
        vals_b,
        color=colors_b,
        edgecolor="black",
        linewidth=0.9,
        width=0.48,
        alpha=0.88,
    )
    for bar in bars_b:
        h = int(bar.get_height())
        pct = 100.0 * h / total if total else 0.0
        ax2.annotate(
            f"{h:,}\n({pct:.1f}%)",
            xy=(bar.get_x() + bar.get_width() / 2, h),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10.5,
            fontweight="bold",
        )
    ax2.set_title(
        f"B. Sample HG00098 Zygosity (Het/Hom Ratio = {het_hom_ratio:.2f})",
        fontsize=13,
        fontweight="bold",
        pad=12,
    )
    ax2.set_ylabel("Variant Count", fontsize=11, fontweight="bold")
    ax2.set_ylim(0, max(vals_b) * 1.24)

    sns.despine(top=True, right=True)
    fig.suptitle(
        "Variant Calling Exploration — Variant Types & Zygosity Profile",
        fontsize=15,
        fontweight="bold",
        y=1.01,
    )
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved variant types & zygosity plot -> {output_png}")


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
        default=ROOT / "results/variant_types.png",
    )
    args = parser.parse_args(argv)
    plot_variant_types(args.input, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
