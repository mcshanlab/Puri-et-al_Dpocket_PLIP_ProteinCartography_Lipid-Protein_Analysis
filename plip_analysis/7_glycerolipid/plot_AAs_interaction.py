import matplotlib
# --- Ensure PDF text stays editable (embed TrueType fonts) ---
matplotlib.rcParams['pdf.fonttype'] = 42   # TrueType
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'   # Force Arial

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

# ----------------------------------
# Load and parse the input text file
# ----------------------------------
input_file = "AA_glycerolipid.txt"

interaction_data = {}
current_category = None

with open(input_file, "r") as f:
    for line in f:
        line = line.strip()

        if line.startswith("###"):
            current_category = line.replace("###", "").strip()
            interaction_data[current_category] = {}
            continue

        m = re.match(r"(\d+)\s+([A-Z]+)", line)
        if m and current_category:
            count = int(m.group(1))
            aa = m.group(2)
            interaction_data[current_category][aa] = count

# ----------------------------------
# Define AA order by biochemical class
# ----------------------------------
groups = {
    "Positively charged": ["ARG", "HIS", "LYS"],
    "Negatively charged": ["ASP", "GLU"],
    "Uncharged": ["SER", "THR", "ASN", "GLN", "CYS"],
    "Aliphatic": ["ALA", "VAL", "ILE", "LEU", "MET", "GLY", "PRO"],
    "Aromatic": ["PHE", "TYR", "TRP"]
}

ordered_aas = [aa for grp in groups.values() for aa in grp]

# ----------------------------------
# Convert parsed text into DataFrame
# ----------------------------------
df = pd.DataFrame(
    0,
    index=ordered_aas,
    columns=interaction_data.keys(),
    dtype=float
)

for category, values in interaction_data.items():
    for aa, count in values.items():
        if aa in df.index:
            df.loc[aa, category] = count

# Normalize each AA row to percentages
df_pct = df.div(df.sum(axis=1), axis=0) * 100

# ----------------------------------
# Plot formatting (small figure)
# ----------------------------------
colors = {
    "Hydrophobic Interactions": "#E88FBC",
    "Hydrogen Bonds": "#7DD3F6",
    "Salt Bridges": "#FFD478",
    "pi-Stacking": "#90C96A",
    "pi-Cation Interactions": "#009152",
    "Halogen Bonds": "#F79420",
    "Metal Complexes": "#C2C0C0"
}

plt.figure(figsize=(4, 2.5))  # Small figure

# Use integer positions for bars
x = np.arange(len(df_pct.index))
bottom = np.zeros(len(df_pct))

for interaction in colors:
    if interaction in df_pct.columns:
        plt.bar(
            x,
            df_pct[interaction],
            bottom=bottom,
            label=interaction,
            color=colors[interaction],
            linewidth=0,
            width=0.65
        )
        bottom += df_pct[interaction].values

# ----------------------------------
# X-axis formatting
# ----------------------------------
plt.xticks(x, df_pct.index, rotation=45, ha="right", fontsize=7)
plt.ylabel("Normalized Count %", fontsize=10)
plt.xlabel("Amino Acids", fontsize=10)
ymin = -5
plt.ylim(ymin, 100)

# ----------------------------------
# Add group labels above bars
# ----------------------------------
ax = plt.gca()
x_positions = {aa: i for i, aa in enumerate(df_pct.index)}

for label, aa_list in groups.items():
    start = x_positions[aa_list[0]]
    end = x_positions[aa_list[-1]]
    mid = (start + end) / 2
    # Place label above tallest bar in the group
    max_height = df_pct.loc[aa_list].sum(axis=1).max()
    ax.text(mid, max_height + 5, label,
            ha='center', va='bottom', fontsize=8, rotation=0)

# ----------------------------------
# Legend outside the plot (right)
# ----------------------------------
plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc='upper left',
    frameon=False,
    fontsize=8
)

# ----------------------------------
# Manual layout adjustment
# ----------------------------------
plt.subplots_adjust(left=0.15, right=0.75, top=0.85, bottom=0.25)

# ----------------------------------
# Save editable PDF
# ----------------------------------
plt.savefig("AA_interactions_glycerolipid.pdf", bbox_inches='tight')
plt.show()
