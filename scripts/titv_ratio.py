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

transitions = {
    ("A", "G"),
    ("G", "A"),
    ("C", "T"),
    ("T", "C")
}

ti = 0
tv = 0

for _, row in df.iterrows():
    ref = row["REF"]
    alt = row["ALT"]

    if len(ref) == 1 and len(alt) == 1:
        if (ref, alt) in transitions:
            ti += 1
        else:
            tv += 1

ratio = ti / tv if tv != 0 else None

print("Transitions:", ti)
print("Transversions:", tv)
print("Ti/Tv ratio:", ratio)