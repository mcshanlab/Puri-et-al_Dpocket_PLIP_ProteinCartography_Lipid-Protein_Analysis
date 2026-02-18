import matplotlib
# Ensure all text in the PDF stays editable (TrueType instead of paths)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Force Arial for all text
matplotlib.rcParams['font.family'] = 'Arial'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# --------------------------------------
# User configuration
# --------------------------------------
base_dir = Path("/Users/amcshan3/Desktop/Manuscripts/PLIP_Dpocket_Lipid_Puri_2025/used-to-make/plip_amino_acid/metal_complex")
output_pdf = "metal_type_barplot.pdf"
counts_output_txt = "metal_type_counts.txt"

# Lipid classes and corresponding subfolders
lipid_dirs = {
    "sterol": "sterol_metalcomplex_stats",
    "Polyketide": "polyketide_metalcomplex_stats",
    "Prenol": "prenol_metalcomplex_stats",
    "saccharo lipid": "saccharolipid_metalcomplex_stats",
    "Spingo lipid": "sphingolipid_metalcomplex_stats",
    "Fatty Acyl": "fattyacyl_metalcomplex_stats",
    "Glycerophopspholipid": "glycerophospholipid_metalcomplex_stats",
    "Glycerolipid": "glycerolipid_metalcomplex_stats"
}

# Custom colors per lipid class
color_map = {
    "sterol": "#F2C9D1",
    "Polyketide": "#CBC7D6",
    "Prenol": "#D0F2F2",
    "saccharo lipid": "#F5F5F5",
    "Spingo lipid": "#F2ECD3",
    "Fatty Acyl": "#F2DDBF",
    "Glycerophopspholipid": "#F5C1CE",
    "Glycerolipid": "#D4EBD1"
}

# --------------------------------------
# Load METAL_TYPE data
# --------------------------------------
all_data = []

for lipid_class, folder in lipid_dirs.items():
    txt_path = base_dir / folder / "METAL_TYPE.txt"

    if not txt_path.exists():
        print(f"WARNING: {txt_path} does not exist, skipping {lipid_class}.")
        continue

    with open(txt_path) as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        print(f"WARNING: {lipid_class} METAL_TYPE.txt is empty. Skipping.")
        continue

    df = pd.DataFrame({"metal_type": lines, "class": lipid_class})
    all_data.append(df)

if not all_data:
    raise ValueError("No valid METAL_TYPE.txt data found.")

merged = pd.concat(all_data, ignore_index=True)

# --------------------------------------
# Count metal types per lipid class
# --------------------------------------
counts = merged.groupby(["class", "metal_type"]).size().unstack(fill_value=0)

# Save raw counts
counts.to_csv(counts_output_txt, sep="\t")
print(f"Saved metal-type counts to: {counts_output_txt}")

# --------------------------------------
# Normalize within each lipid class
# --------------------------------------
normalized = counts.div(counts.max(axis=1).replace(0, np.nan), axis=0).fillna(0)

# --------------------------------------
# Plot bar chart using normalized counts
# --------------------------------------
plt.figure(figsize=(10, 5))

x = np.arange(len(normalized.index))
width = 0.75 / len(normalized.columns)   # bar width distribution

for i, metal in enumerate(normalized.columns):
    plt.bar(
        x + i * width,
        normalized[metal],
        width=width,
        label=metal
    )

plt.xticks(
    x + width * (len(normalized.columns)-1) / 2,
    normalized.index,
    rotation=40,
    ha="right"
)

plt.ylabel("Normalized count (relative to largest in class)", fontsize=11)
plt.title("Normalized Metal-Type Frequencies per Lipid Class (max = 1.0)")

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.ylim(0, 1.05)

# --------------------------------------
# Legend outside the plot
# --------------------------------------
plt.legend(
    title="Metal type",
    fontsize=8,
    loc="center left",
    bbox_to_anchor=(1, 0.5)
)

plt.tight_layout(rect=[0, 0, 0.85, 1])  # leave space for legend
plt.savefig(output_pdf, dpi=300)
plt.close()

print(f"Saved normalized figure to: {output_pdf}")
