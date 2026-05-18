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

plt.bar(["Transition", "Transversion"], [ti, tv])
plt.title("Ti/Tv distribution")
plt.savefig("results/titv_plot.png")