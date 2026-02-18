#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# USER SETTINGS
# ============================================================
INPUT_CSV = "ss_by_LeidenCluster.csv"
OUTPUT_PDF = "ss_by_LeidenCluster.pdf"

# ============================================================
# LOAD DATA
# ============================================================
df = pd.read_csv(INPUT_CSV)
df = df.sort_values("LeidenCluster")

labels = df["LeidenCluster"].values
x = np.arange(len(labels))

helix_mean = df["helix_mean"].values
sheet_mean = df["sheet_mean"].values
loop_mean = df["loop_mean"].values

# ============================================================
# PLOT STACKED BARS
# ============================================================
fig, ax = plt.subplots(figsize=(12, 6))

ax.bar(x, helix_mean, label="Helix", color="#1f77b4")
ax.bar(x, sheet_mean, bottom=helix_mean, label="Sheet", color="#ff7f0e")
ax.bar(x, loop_mean, bottom=helix_mean + sheet_mean, label="Loop", color="#2ca02c")

# ------------------------------------------------------------
# AXES & LABELS
# ------------------------------------------------------------
ax.set_xticks(x)
ax.set_xticklabels(
    labels,
    rotation=45,
    ha="right",
    rotation_mode="anchor",
    fontsize=9
)

# Give extra room on the left for rotated labels
ax.set_xlim(-0.5, len(x) - 0.5)

# MANUAL layout control (do NOT use tight_layout)
plt.subplots_adjust(left=0.18, right=0.82, bottom=0.25)

ax.set_ylabel("Percentage of residues near ligand (%)")
ax.set_xlabel("LeidenCluster")
ax.set_title("Secondary Structure Composition near Ligand by Leiden Cluster")

# ------------------------------------------------------------
# LEGEND OUTSIDE
# ------------------------------------------------------------
ax.legend(
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    frameon=False
)

plt.tight_layout(rect=[0, 0, 0.85, 1])

# ============================================================
# SAVE
# ============================================================
plt.savefig(OUTPUT_PDF, format="pdf")
print(f"Plot saved as editable PDF: {OUTPUT_PDF}")
plt.show()
