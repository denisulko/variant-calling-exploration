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

def classify(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    return "INDEL"

df["TYPE"] = df.apply(lambda row: classify(row["REF"], row["ALT"]), axis=1)

counts = df["TYPE"].value_counts()

counts.plot(kind="bar")
plt.title("Variant types")
plt.savefig("results/variant_types.png")