from pathlib import Path

file_path = Path("data/raw/sample.vcf")

variants = []

with open(file_path, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        variants.append(line.strip().split("\t"))

print("Total variants:", len(variants))
print("First variant:", variants[0][:5])