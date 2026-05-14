import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv("polyketides_plip_interaction_counts.tsv", sep="\t")

# -----------------------------
# Define NEW plotting order
# -----------------------------
order = [
    "hydrophobic",
    "hbond",
    "pistacking",
    "pication",
    "halogen",
    "saltbridge",
    "metal"
]

labels = [
    "Hydrophobic",
    "Hydrogen bond",
    "π-stacking",
    "π-cation",
    "Halogen bond",
    "Salt bridge",
    "Metal complex"
]

# -----------------------------
# Sum across all structures
# -----------------------------
raw_values = np.array([df[col].sum() for col in order])

# -----------------------------
# NORMALIZATION (max = 100)
# -----------------------------
max_val = raw_values.max()
values = (raw_values / max_val) * 100

# -----------------------------
# Colors
# -----------------------------
colors = [
    "#E88FBC",  # hydrophobic
    "#7DD3F6",  # hydrogen bond
    "#90C96A",  # π-stacking
    "#009152",  # π-cation
    "#F79420",  # halogen bond
    "#FFD478",  # salt bridge
    "#C2C0C0"   # metal complex
]

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(10, 6))

bars = plt.bar(labels, values, color=colors)

plt.xticks(rotation=30, ha="right")
plt.ylabel("Normalized number of interactions")
plt.title("Polyketides–Protein Interaction Profiles (Normalized)")

# -----------------------------
# Legend
# -----------------------------
legend_handles = [
    plt.Rectangle((0, 0), 1, 1, color=c)
    for c in colors
]

plt.legend(
    legend_handles,
    labels,
    title="Interaction Types",
    bbox_to_anchor=(1.05, 1),
    loc="upper left"
)

plt.tight_layout()

# -----------------------------
# Save PDF
# -----------------------------
plt.savefig(
    "polyketides_interaction_profiles_normalized.pdf",
    format="pdf",
    bbox_inches="tight"
)

plt.show()
