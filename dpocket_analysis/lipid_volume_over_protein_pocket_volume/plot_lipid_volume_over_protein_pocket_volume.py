import matplotlib
# Ensure all text in the PDF stays editable (TrueType instead of paths)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Force Arial for all text
matplotlib.rcParams['font.family'] = 'Arial'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --------------------------------------
# User configuration
# --------------------------------------
excel_file = "../../Source_Data_dpocket.xlsx"

target_sheets = [
    "Sterol",
    "Polyketide",
    "Prenol",
    "Saccharolipid",
    "Sphingolipid",
    "Fatty Acyl",
    "Glycerophospholipid",
    "Glycerolipid"
]

lig_col_keyword = "lig_vol"
pocket_col_keyword = "pock_vol"

output_pdf = "lipid_to_pocket_volume_ratio_violin.pdf"
median_output_txt = "lipid_to_pocket_volume_ratio_medians.txt"

# Custom colors per lipid class
color_map = {
    "Sterol": "#F2C9D1",
    "Polyketide": "#CBC7D6",
    "Prenol": "#D0F2F2",
    "Saccharolipid": "#F5F5F5",
    "Sphingolipid": "#F2ECD3",
    "Fatty Acyl": "#F2DDBF",
    "Glycerophospholipid": "#F5C1CE",
    "Glycerolipid": "#D4EBD1"
}

# --------------------------------------
# Load data
# --------------------------------------
all_data = []

for sheet in target_sheets:
    df = pd.read_excel(excel_file, sheet_name=sheet)

    lig_col = None
    pock_col = None

    # auto-find columns
    for col in df.columns:
        if lig_col_keyword in col.lower():
            lig_col = col
        if pocket_col_keyword in col.lower():
            pock_col = col

    if lig_col is None or pock_col is None:
        print(f"WARNING: required columns not found in sheet {sheet}, skipping...")
        continue

    tmp = df[[lig_col, pock_col]].copy()
    tmp[lig_col] = pd.to_numeric(tmp[lig_col], errors='coerce')
    tmp[pock_col] = pd.to_numeric(tmp[pock_col], errors='coerce')

    # remove zeros to avoid divide-by-zero
    tmp = tmp.dropna(subset=[lig_col, pock_col])
    tmp = tmp[tmp[pock_col] > 0]

    if tmp.empty:
        print(f"WARNING: sheet {sheet} has no valid ratio data, skipping...")
        continue

    # Calculate ratio
    tmp["ratio"] = tmp[lig_col] / tmp[pock_col]
    tmp["class"] = sheet

    all_data.append(tmp[["ratio", "class"]])

if not all_data:
    raise ValueError("No valid ratio data found in any sheet.")

merged = pd.concat(all_data, ignore_index=True)
merged["class"] = pd.Categorical(merged["class"], categories=target_sheets, ordered=True)

# --------------------------------------
# Prepare data for violin plot
# --------------------------------------
data_list = [merged[merged["class"] == c]["ratio"] for c in target_sheets]
sample_sizes = [len(vals) for vals in data_list]

# --------------------------------------
# Calculate medians and save to text file
# --------------------------------------
medians = [vals.median() for vals in data_list]

with open(median_output_txt, "w") as f:
    f.write("Lipid class\tMedian (Lipid Volume / Pocket Volume)\n")
    for cls, med in zip(target_sheets, medians):
        f.write(f"{cls}\t{med:.3f}\n")

print(f"Saved median values to: {median_output_txt}")

# --------------------------------------
# Plot
# --------------------------------------
plt.figure(figsize=(4.5, 4))

vp = plt.violinplot(
    data_list,
    showmeans=False,
    showextrema=False,
    showmedians=True
)

# Apply colors
for body, cls in zip(vp["bodies"], target_sheets):
    body.set_facecolor(color_map[cls])
    body.set_edgecolor("black")
    body.set_linewidth(1)
    body.set_alpha(0.95)

vp["cmedians"].set_color("black")
vp["cmedians"].set_linewidth(1.5)

# Axes formatting
plt.xticks(range(1, len(target_sheets) + 1), target_sheets, rotation=40, ha="right")
plt.ylabel("Lipid Volume / Protein Pocket Volume", fontsize=11)

ax = plt.gca()
ax.set_ylim(0, 3)  # <-- force y-axis value if you want
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_linewidth(1)

# Add sample sizes
ymax = merged["ratio"].max()

for i, n in enumerate(sample_sizes, start=1):
    plt.text(
        i,
        ymax * 1.05,
        f"n = {n}",
        ha="center",
        va="bottom",
        fontsize=9,
        rotation=45
    )

plt.tight_layout()
plt.savefig(output_pdf, dpi=300)
plt.close()

print(f"Saved figure to: {output_pdf}")
