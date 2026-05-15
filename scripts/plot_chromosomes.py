from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

file_path = Path("data/raw/sample.vcf")

records = []

with open(file_path, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        records.append(parts[:5])

df = pd.DataFrame(records, columns=["CHROM", "POS", "ID", "REF", "ALT"])

counts = df["CHROM"].value_counts()

counts.plot(kind="bar")
plt.xlabel("Chromosome")
plt.ylabel("Variant count")
plt.title("Variants per chromosome")

plt.tight_layout()
plt.savefig("results/chromosome_distribution.png")