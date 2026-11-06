"""Plot Transition/Transversion (Ti/Tv) ratio and 12-class substitution spectrum."""

from __future__ import annotations

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

from read_vcf import ROOT, TRANSITION_PAIRS
from titv_ratio import compute_titv_metrics


def plot_titv_distribution(vcf_path: Path, output_png: Path) -> None:
    """Render two-panel Ti/Tv summary and 12-substitution spectrum."""
    metrics = compute_titv_metrics(vcf_path)
    ti = int(metrics["transitions"])
    tv = int(metrics["transversions"])
    ratio = float(metrics["titv_ratio"])
    pair_counts = metrics["pair_counts"]
    total_snps = ti + tv

    sns.set_theme(style="whitegrid", context="talk")
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(15, 6.2), dpi=300, gridspec_kw={"width_ratios": [1, 1.55]}
    )

    # Panel A: Transitions vs Transversions
    bars1 = ax1.bar(
        ["Transitions (Ti)\n(Purine↔Purine / Pyr↔Pyr)", "Transversions (Tv)\n(Purine↔Pyrimidine)"],
        [ti, tv],
        color=["#1f77b4", "#e67e22"],
        edgecolor="black",
        linewidth=1.0,
        width=0.52,
        alpha=0.88,
    )
    for bar in bars1:
        h = int(bar.get_height())
        pct = 100.0 * h / total_snps if total_snps else 0.0
        ax1.annotate(
            f"{h:,}\n({pct:.1f}%)",
            xy=(bar.get_x() + bar.get_width() / 2, h),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10.5,
            fontweight="bold",
        )
    ax1.set_title(
        f"A. Ti vs. Tv Count (Observed Ti/Tv = {ratio:.2f})",
        fontsize=13,
        fontweight="bold",
        pad=12,
    )
    ax1.set_ylabel("SNP Count", fontsize=11, fontweight="bold")
    ax1.set_ylim(0, max(ti, tv) * 1.25)

    # Panel B: 12-class single-nucleotide substitution spectrum
    ordered_pairs = sorted(pair_counts.keys(), key=lambda p: pair_counts[p], reverse=True)
    pair_labels = [f"{r}>{a}" for r, a in ordered_pairs]
    pair_vals = [pair_counts[p] for p in ordered_pairs]
    pair_colors = ["#1f77b4" if p in TRANSITION_PAIRS else "#e67e22" for p in ordered_pairs]

    bars2 = ax2.bar(
        pair_labels,
        pair_vals,
        color=pair_colors,
        edgecolor="black",
        linewidth=0.8,
        alpha=0.88,
    )
    for bar in bars2:
        h = int(bar.get_height())
        ax2.annotate(
            f"{h:,}",
            xy=(bar.get_x() + bar.get_width() / 2, h),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8.5,
            fontweight="bold",
        )
    ax2.set_title(
        "B. 12-Class Single-Nucleotide Substitution Spectrum (CpG Enrichment)",
        fontsize=13,
        fontweight="bold",
        pad=12,
    )
    ax2.set_xlabel("Base Substitution (REF > ALT)", fontsize=11, fontweight="bold")
    ax2.set_ylabel("SNP Count", fontsize=11, fontweight="bold")
    ax2.set_ylim(0, max(pair_vals) * 1.2)

    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="#1f77b4", edgecolor="black", label="Transition (4 pathways: A↔G, C↔T)"),
        Patch(facecolor="#e67e22", edgecolor="black", label="Transversion (8 pathways)"),
    ]
    ax2.legend(handles=legend_handles, frameon=True, facecolor="white", fontsize=9.5, loc="upper right")

    sns.despine(top=True, right=True)
    fig.suptitle(
        f"Variant Calling QC — Transition / Transversion (Ti/Tv = {ratio:.2f}) Analysis",
        fontsize=15,
        fontweight="bold",
        y=1.01,
    )
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved Ti/Tv distribution figure -> {output_png}")


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
        default=ROOT / "results/titv_plot.png",
    )
    args = parser.parse_args(argv)
    plot_titv_distribution(args.input, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
