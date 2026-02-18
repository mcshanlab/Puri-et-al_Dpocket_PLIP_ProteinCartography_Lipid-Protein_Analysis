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
input_dir = "/Users/amcshan3/Desktop/Manuscripts/PLIP_Dpocket_Lipid_Puri_2025/used-to-make/plip_amino_acid/hydrophobic_distances"
output_pdf = "hydrophobic_distance_violin.pdf"
median_output_txt = "hydrophobic_distance_median_values.txt"

# Lipid classes and corresponding input files
lipid_files = {
    "sterol": "hydrophobic_distances_sterol.txt",
    "Polyketide": "hydrophobic_distances_polyketide.txt",
    "Prenol": "hydrophobic_distances_prenol.txt",
    "saccharo lipid": "hydrophobic_distances_saccharolipid.txt",
    "Spingo lipid": "hydrophobic_distances_sphingolipid.txt",
    "Fatty Acyl": "hydrophobic_distances_fattyacyl.txt",
    "Glycerophopspholipid": "hydrophobic_distances_glycerophospholipid.txt",
    "Glycerolipid": "hydrophobic_distances_glycerolipid.txt"
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
# Load data
# --------------------------------------
all_data = []

for lipid_class, filename in lipid_files.items():
    path = Path(input_dir) / filename
    if not path.exists():
        print(f"WARNING: file {filename} not found, skipping {lipid_class}")
        continue

    # Skip header line if present
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip() and not line.lower().startswith("hydrophobic")]

    if not lines:
        print(f"WARNING: {filename} has no data, skipping {lipid_class}")
        continue

    distances = pd.to_numeric(pd.Series(lines), errors='coerce').dropna()
    if distances.empty:
        print(f"WARNING: {filename} has no numeric values, skipping {lipid_class}")
        continue

    df = pd.DataFrame({
        "distance": distances,
        "class": lipid_class
    })
    all_data.append(df)

if not all_data:
    raise ValueError("No valid distance data found in any file.")

merged = pd.concat(all_data, ignore_index=True)

# --------------------------------------
# Only include classes that have data
# --------------------------------------
target_classes = merged["class"].unique().tolist()
merged["class"] = pd.Categorical(merged["class"], categories=target_classes, ordered=True)

# Prepare data for violin plot
data_list = [merged[merged["class"] == c]["distance"] for c in target_classes]

# Sample sizes
sample_sizes = [len(vals) for vals in data_list]

# --------------------------------------
# Calculate medians and save to text file
# --------------------------------------
medians = [vals.median() for vals in data_list]

with open(median_output_txt, "w") as f:
    f.write("Lipid class\tMedian hydrophobic distance (Å)\n")
    for cls, med in zip(target_classes, medians):
        f.write(f"{cls}\t{med:.3f}\n")

print(f"Saved median values to: {median_output_txt}")

# --------------------------------------
# Plot
# --------------------------------------
plt.figure(figsize=(6, 4))

vp = plt.violinplot(
    data_list,
    showmeans=False,
    showextrema=False,
    showmedians=True
)

# Apply unique colors per lipid class
for body, cls in zip(vp["bodies"], target_classes):
    body.set_facecolor(color_map.get(cls, "#CCCCCC"))
    body.set_edgecolor("black")
    body.set_linewidth(1)
    body.set_alpha(0.95)

# Median line style
vp["cmedians"].set_color("black")
vp["cmedians"].set_linewidth(1.5)

# Axes formatting
plt.xticks(range(1, len(target_classes) + 1), target_classes, rotation=40, ha="right")
plt.ylabel("Hydrophobic interaction distance (Å)", fontsize=11)

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_linewidth(1)

# ------------------------------------------------------
# NEW: Add median ± SD labels above each violin
# ------------------------------------------------------
ymax = merged["distance"].max()
for i, vals in enumerate(data_list, start=1):
    med = vals.median()
    sd = vals.std()
    plt.text(
        i,
        ymax * 1.05,
        f"{med:.2f} ± {sd:.2f} Å",
        ha="center",
        va="bottom",
        fontsize=9,
        rotation=45
    )

plt.tight_layout()

# Save as PDF with editable text
plt.savefig(output_pdf, dpi=300)
plt.close()

print(f"Saved figure to: {output_pdf}")
