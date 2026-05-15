from pathlib import Path
import pandas as pd

file_path = Path("data/raw/sample.vcf")

records = []

with open(file_path, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        records.append(parts[:5])

df = pd.DataFrame(records, columns=["CHROM", "POS", "ID", "REF", "ALT"])

print(df.head())
print()
print(df["CHROM"].value_counts())

df.to_csv("results/vcf_summary.csv", index=False)