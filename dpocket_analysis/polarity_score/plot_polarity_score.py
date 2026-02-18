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

# NEW: column for proportion of polar atoms
prop_col_keyword = "polarity_score"
output_pdf = "polarity_score_violin.pdf"
median_output_txt = "polarity_score_median_values.txt"

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

    # auto-find polarity_score column
    prop_col = None
    for col in df.columns:
        if prop_col_keyword in col.lower():
            prop_col = col
            break

    if prop_col is None:
        print(f"WARNING: no '{prop_col_keyword}' column found in sheet {sheet}, skipping...")
        continue

    tmp = df[[prop_col]].copy()
    tmp[prop_col] = pd.to_numeric(tmp[prop_col], errors='coerce')
    tmp["class"] = sheet
    tmp = tmp.dropna(subset=[prop_col])

    if tmp.empty:
        print(f"WARNING: sheet {sheet} has no numeric values for '{prop_col_keyword}', skipping...")
        continue

    all_data.append(tmp[[prop_col, "class"]])

if not all_data:
    raise ValueError(f"No valid '{prop_col_keyword}' data found in any sheet.")

merged = pd.concat(all_data, ignore_index=True)

# order categories for consistent plotting
merged["class"] = pd.Categorical(merged["class"], categories=target_sheets, ordered=True)

# --------------------------------------
# Prepare data for violin plot
# --------------------------------------
data_list = [merged[merged["class"] == c][prop_col] for c in target_sheets]

# sample sizes
sample_sizes = [len(vals) for vals in data_list]

# --------------------------------------
# Calculate medians and save to text file
# --------------------------------------
medians = [vals.median() for vals in data_list]

with open(median_output_txt, "w") as f:
    f.write("Polarity score\n")
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

# --------------------------------------
# Apply unique colors per lipid class
# --------------------------------------
for body, cls in zip(vp["bodies"], target_sheets):
    body.set_facecolor(color_map[cls])
    body.set_edgecolor("black")
    body.set_linewidth(1)
    body.set_alpha(0.95)

# Median line style
vp["cmedians"].set_color("black")
vp["cmedians"].set_linewidth(1.5)

# Axes formatting
plt.xticks(range(1, len(target_sheets) + 1), target_sheets, rotation=40, ha="right")
plt.ylabel("Polarity score", fontsize=11)

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_linewidth(1)

# --------------------------------------
# Add n = X labels above each violin
# --------------------------------------
ymax = merged[prop_col].max()

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
