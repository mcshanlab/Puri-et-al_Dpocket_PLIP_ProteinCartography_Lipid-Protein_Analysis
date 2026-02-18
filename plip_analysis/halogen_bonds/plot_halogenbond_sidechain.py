import matplotlib
# Ensure all text in the PDF stays editable (TrueType instead of paths)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# --------------------------------------
# User configuration
# --------------------------------------
base_dir = Path("/Users/amcshan3/Desktop/Manuscripts/PLIP_Dpocket_Lipid_Puri_2025/used-to-make/plip_amino_acid/halogen_bonds")
output_pdf = "sidechain_counts_barplot.pdf"
output_txt = "sidechain_counts_normalized.txt"

# Lipid classes and corresponding subfolders
lipid_dirs = {
    "sterol": "sterol_halogenbond_stats",
    "Polyketide": "polyketide_halogenbond_stats",
    "Prenol": "prenol_halogenbond_stats",
    "saccharo lipid": "saccharolipid_halogenbond_stats",
    "Spingo lipid": "sphingolipid_halogenbond_stats",
    "Fatty Acyl": "fattyacyl_halogenbond_stats",
    "Glycerophopspholipid": "glycerophospholipid_halogenbond_stats",
    "Glycerolipid": "glycerolipid_halogenbond_stats"
}

# Colors for True / False
color_map = {
    "True": "#4C72B0",   # side-chain
    "False": "#55A868"   # backbone
}

# --------------------------------------
# Load and count data
# --------------------------------------
all_data = []

for lipid_class, folder in lipid_dirs.items():
    txt_path = base_dir / folder / "sidechain.txt"

    if not txt_path.exists():
        print(f"WARNING: {txt_path} does not exist, skipping {lipid_class}.")
        continue

    with open(txt_path) as f:
        lines = [line.strip() for line in f if line.strip()]

    # Convert to boolean
    bools = pd.Series([line.lower() == "true" for line in lines])

    if bools.empty:
        print(f"WARNING: {lipid_class} has no valid True/False entries. Skipping.")
        continue

    counts = bools.value_counts().reindex([True, False], fill_value=0)
    counts.name = lipid_class
    all_data.append(counts)

if not all_data:
    raise ValueError("No valid sidechain.txt data found.")

df_counts = pd.DataFrame(all_data)

# --------------------------------------
# Normalize within each lipid class
# --------------------------------------
df_counts_norm = df_counts.div(df_counts.max(axis=1), axis=0)
df_counts_norm.index.name = "class"
df_counts_norm.reset_index(inplace=True)

# Save normalized counts
df_counts_norm.to_csv(output_txt, sep="\t", index=False)
print(f"Saved normalized counts to: {output_txt}")

# --------------------------------------
# Plot bar chart with legend outside
# --------------------------------------
fig, ax = plt.subplots(figsize=(10, 5))

x = df_counts_norm["class"]
bar_width = 0.35
indices = range(len(x))

bars_true = ax.bar(indices, df_counts_norm[True], bar_width, label="Side-chain (True)", color=color_map["True"])
bars_false = ax.bar([i + bar_width for i in indices], df_counts_norm[False], bar_width, label="Backbone (False)", color=color_map["False"])

ax.set_xticks([i + bar_width / 2 for i in indices])
ax.set_xticklabels(x, rotation=45, ha="right")
ax.set_ylabel("Normalized Halogen Counts")
ax.set_title("Normalized Counts of Side-chain vs Backbone Halogen Bonds per Lipid Class")

# Move legend outside
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

plt.tight_layout(rect=[0, 0, 0.85, 1])  # leave space on the right for the legend
plt.savefig(output_pdf, dpi=300)
plt.close()
print(f"Saved figure to: {output_pdf}")
