import matplotlib
# Ensure PDF text stays editable (TrueType)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# --------------------------------------
# User configuration
# --------------------------------------
base_dir = Path("/Users/amcshan3/Desktop/Manuscripts/PLIP_Dpocket_Lipid_Puri_2025/used-to-make/plip_amino_acid/pi-stacking")
output_pdf = "pistacking_angle_violin.pdf"
median_output_txt = "pistacking_angle_stats.txt"

# Lipid classes and corresponding subfolders
lipid_dirs = {
    "sterol": "sterol_pistacking_stats",
    "Polyketide": "polyketide_pistacking_stats",
    "Prenol": "prenol_pistacking_stats",
    "saccharo lipid": "saccharolipid_pistacking_stats",
    "Spingo lipid": "sphingolipid_pistacking_stats",
    "Fatty Acyl": "fattyacyl_pistacking_stats",
    "Glycerophopspholipid": "glycerophospholipid_pistacking_stats",
    "Glycerolipid": "glycerolipid_pistacking_stats"
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
# Load ANGLE data
# --------------------------------------
all_data = []

for lipid_class, folder in lipid_dirs.items():
    txt_path = base_dir / folder / "ANGLE.txt"
    if not txt_path.exists():
        print(f"WARNING: {txt_path} does not exist, skipping {lipid_class}.")
        continue

    with open(txt_path) as f:
        lines = [line.strip() for line in f if line.strip()]

    angles = pd.to_numeric(pd.Series(lines), errors="coerce").dropna()
    if angles.empty:
        print(f"WARNING: {lipid_class} has no numeric ANGLE values. Skipping.")
        continue

    df = pd.DataFrame({"angle": angles, "class": lipid_class})
    all_data.append(df)

if not all_data:
    raise ValueError("No valid ANGLE.txt data found.")

merged = pd.concat(all_data, ignore_index=True)

# --------------------------------------
# Prepare for plotting
# --------------------------------------
target_classes = merged["class"].unique().tolist()
merged["class"] = pd.Categorical(merged["class"], categories=target_classes, ordered=True)

data_list = [merged[merged["class"] == c]["angle"] for c in target_classes]

# --------------------------------------
# Compute median ± SD
# --------------------------------------
medians = [vals.median() for vals in data_list]
sds = [vals.std() for vals in data_list]

with open(median_output_txt, "w") as f:
    f.write("Lipid class\tMedian (°)\tSD (°)\tN\n")
    for cls, med, sd, vals in zip(target_classes, medians, sds, data_list):
        f.write(f"{cls}\t{med:.2f}\t{sd:.2f}\t{len(vals)}\n")

print(f"Saved stats to: {median_output_txt}")

# --------------------------------------
# Plot violin
# --------------------------------------
plt.figure(figsize=(10, 5))

vp = plt.violinplot(
    data_list,
    showmeans=False,
    showextrema=False,
    showmedians=True
)

# Apply colors per lipid class
for body, cls in zip(vp["bodies"], target_classes):
    body.set_facecolor(color_map.get(cls, "#CCCCCC"))
    body.set_edgecolor("black")
    body.set_linewidth(1)
    body.set_alpha(0.95)

vp["cmedians"].set_color("black")
vp["cmedians"].set_linewidth(1.5)

# Axes formatting
plt.xticks(range(1, len(target_classes) + 1), target_classes, rotation=40, ha="right")
plt.ylabel("π-Stacking Angle (°)", fontsize=11)

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Median ± SD labels above violins
ymax = max([vals.max() for vals in data_list])
for i, (med, sd) in enumerate(zip(medians, sds), start=1):
    plt.text(
        i,
        ymax * 1.05,
        f"{med:.2f} ± {sd:.2f}°",
        ha="center",
        va="bottom",
        fontsize=9,
        rotation=45
    )

plt.tight_layout()
plt.savefig(output_pdf, dpi=300)
plt.close()

print(f"Saved figure to: {output_pdf}")
