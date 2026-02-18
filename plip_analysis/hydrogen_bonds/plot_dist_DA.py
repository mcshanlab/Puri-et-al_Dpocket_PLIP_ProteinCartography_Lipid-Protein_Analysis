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
base_dir = Path("/Users/amcshan3/Desktop/Manuscripts/PLIP_Dpocket_Lipid_Puri_2025/used-to-make/plip_amino_acid/hydrogen_bonds")
output_pdf = "donor_acceptor_distance_violin.pdf"
median_output_txt = "donor_acceptor_distance_medians.txt"

# Lipid classes and corresponding subfolders
lipid_dirs = {
    "sterol": "sterol_hbond_stats",
    "Polyketide": "polyketide_hbond_stats",
    "Prenol": "prenol_hbond_stats",
    "saccharo lipid": "saccharolipid_hbond_stats",
    "Spingo lipid": "sphingolipid_hbond_stats",
    "Fatty Acyl": "fattyacyl_hbond_stats",
    "Glycerophopspholipid": "glycerophospholipid_hbond_stats",
    "Glycerolipid": "glycerolipid_hbond_stats"
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

for lipid_class, folder in lipid_dirs.items():
    txt_path = base_dir / folder / "dist_DA.txt"

    if not txt_path.exists():
        print(f"WARNING: {txt_path} does not exist, skipping {lipid_class}.")
        continue

    with open(txt_path) as f:
        lines = [line.strip() for line in f if line.strip()]

    distances = pd.to_numeric(pd.Series(lines), errors="coerce").dropna()

    if distances.empty:
        print(f"WARNING: {lipid_class} has no numeric donor distances. Skipping.")
        continue

    df = pd.DataFrame({"distance": distances, "class": lipid_class})
    all_data.append(df)

if not all_data:
    raise ValueError("No valid dist_DA.txt data found.")

merged = pd.concat(all_data, ignore_index=True)

# --------------------------------------
# Prepare for plotting
# --------------------------------------
target_classes = merged["class"].unique().tolist()
merged["class"] = pd.Categorical(merged["class"], categories=target_classes, ordered=True)

data_list = [merged[merged["class"] == c]["distance"] for c in target_classes]

# --------------------------------------
# Save medians
# --------------------------------------
medians = [vals.median() for vals in data_list]
sds = [vals.std() for vals in data_list]

with open(median_output_txt, "w") as f:
    f.write("Lipid class\tMedian (Å)\tSD (Å)\n")
    for cls, med, sd in zip(target_classes, medians, sds):
        f.write(f"{cls}\t{med:.3f}\t{sd:.3f}\n")

print(f"Saved median values to: {median_output_txt}")

# --------------------------------------
# Plot violin
# --------------------------------------
plt.figure(figsize=(6, 4))

vp = plt.violinplot(
    data_list,
    showmeans=False,
    showextrema=False,
    showmedians=True
)

# Apply colors
for body, cls in zip(vp["bodies"], target_classes):
    body.set_facecolor(color_map.get(cls, "#CCCCCC"))
    body.set_edgecolor("black")
    body.set_linewidth(1)
    body.set_alpha(0.95)

vp["cmedians"].set_color("black")
vp["cmedians"].set_linewidth(1.5)

# Axis formatting
plt.xticks(range(1, len(target_classes) + 1), target_classes, rotation=40, ha="right")
plt.ylabel("Donor–Acceptor Distance (Å)", fontsize=11)

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Labels above violins
ymax = merged["distance"].max()
for i, (med, sd) in enumerate(zip(medians, sds), start=1):
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
plt.savefig(output_pdf, dpi=300)
plt.close()

print(f"Saved figure to: {output_pdf}")
