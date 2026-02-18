import matplotlib
# Ensure editable text in the PDF (TrueType)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
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

pock_col_keyword = "pock"   # used to find the pock_vol column
output_pdf = "pocket_volume_violin.pdf"
median_output_txt = "pocket_volumes_median_values.txt"

# Custom colors per lipid class (same palette as before)
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

    # auto-find pock_vol column
    pock_col = None
    for col in df.columns:
        if pock_col_keyword in col.lower():
            pock_col = col
            break

    if pock_col is None:
        print(f"WARNING: no pocket-volume column found in sheet {sheet}")
        continue

    tmp = df[[pock_col]].copy()
    tmp.columns = ["pock_vol"]
    tmp["pock_vol"] = pd.to_numeric(tmp["pock_vol"], errors="coerce")
    tmp = tmp.dropna(subset=["pock_vol"])
    tmp["class"] = sheet

    all_data.append(tmp)

merged = pd.concat(all_data, ignore_index=True)

# Maintain ordering
merged["class"] = pd.Categorical(merged["class"], categories=target_sheets, ordered=True)

# --------------------------------------
# Prepare data for violin plot
# --------------------------------------
data_list = [merged[merged["class"] == c]["pock_vol"] for c in target_sheets]

# Sample sizes
sample_sizes = [len(vals) for vals in data_list]

# --------------------------------------
# Calculate medians and save
# --------------------------------------
medians = [vals.median() for vals in data_list]

with open(median_output_txt, "w") as f:
    f.write("Lipid class\tMedian pocket volume (Å³)\n")
    for cls, med in zip(target_sheets, medians):
        f.write(f"{cls}\t{med:.2f}\n")

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

# Color each violin individually
for body, cls in zip(vp["bodies"], target_sheets):
    body.set_facecolor(color_map[cls])
    body.set_edgecolor("black")
    body.set_linewidth(1)
    body.set_alpha(0.95)

# Median line
vp["cmedians"].set_color("black")
vp["cmedians"].set_linewidth(1.5)

# Axes formatting
plt.xticks(range(1, len(target_sheets) + 1), target_sheets, rotation=40, ha="right")
plt.ylabel("Pocket volume (Å³)", fontsize=11)

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_linewidth(1)

# --------------------------------------
# Add n = X labels
# --------------------------------------
ymax = merged["pock_vol"].max()

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

# --------------------------------------
# Save as PDF with editable text
# --------------------------------------
plt.savefig(output_pdf, dpi=300)
plt.close()

print(f"Saved figure to: {output_pdf}")
