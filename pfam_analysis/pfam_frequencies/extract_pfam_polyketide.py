#!/usr/bin/env python3
"""
Extract PFAM summary (one PFAM per BioDolphinID), count unique BioDolphinIDs per PFAM,
compute frequency relative to total unique BioDolphinIDs, and plot pie chart for PFAMs >1%.

Outputs:
- polyketide_pfams_summary.xlsx (BioDolphinID, protein_Pfam_ID, count, frequency_percent)
- polyketide_pfam_pie_chart.pdf (editable text, 50% transparent wedges)
"""

import ast
import re
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import patheffects

# ---------------------------
# Ensure editable text in PDFs
# ---------------------------
matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'

# ---------------------------
# Input file
# ---------------------------
excel_path = "/Volumes/8TB_McShan_Drive/nikki/project_1/duplicates_deleted/polyketide.xlsx"
df = pd.read_excel(excel_path)

# ---------------------------
# Validate required columns
# ---------------------------
required_cols = ["BioDolphinID", "protein_Pfam_ID"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"ERROR: Missing required columns: {missing}")

# ---------------------------
# Helper to parse PFAM cell into a list of PFAM IDs
# ---------------------------
def parse_pfam_cell(x):
    if pd.isna(x):
        return []
    if isinstance(x, (list, tuple, set)):
        return [str(s).strip() for s in x if str(s).strip()]
    s = str(x).strip()
    if s.startswith('[') and s.endswith(']'):
        try:
            parsed = ast.literal_eval(s)
            if isinstance(parsed, (list, tuple, set)):
                return [str(ii).strip() for ii in parsed if str(ii).strip()]
        except Exception:
            pass
    if ',' in s or ';' in s:
        tokens = re.split(r'[;,]\s*', s)
        return [t.strip() for t in tokens if t.strip()]
    if s:
        return [s]
    return []

# ---------------------------
# Expand rows: one PFAM per row per BioDolphinID
# ---------------------------
rows = []
for _, row in df.iterrows():
    bd = row["BioDolphinID"]
    pfam_cell = row["protein_Pfam_ID"]
    pfams = parse_pfam_cell(pfam_cell)
    for p in pfams:
        rows.append({"BioDolphinID": bd, "protein_Pfam_ID": p})

expanded = pd.DataFrame(rows)
if expanded.empty:
    raise ValueError("No PFAM entries found after parsing. Check protein_Pfam_ID formatting.")

# ---------------------------
# Count unique BioDolphinIDs per PFAM
# ---------------------------
pfam_counts = (
    expanded.groupby("protein_Pfam_ID")["BioDolphinID"]
    .nunique()
    .reset_index(name="count")
)

total_unique_biodolphin = expanded["BioDolphinID"].nunique()
pfam_counts["frequency_percent"] = pfam_counts["count"] / total_unique_biodolphin * 100

# ---------------------------
# Merge counts onto unique BioDolphin–PFAM pairs
# ---------------------------
unique_pairs = expanded.drop_duplicates(subset=["BioDolphinID", "protein_Pfam_ID"])
result = unique_pairs.merge(pfam_counts, on="protein_Pfam_ID", how="left")

# Sort so highest-frequency PFAMs first
result = result.sort_values(by=["frequency_percent", "protein_Pfam_ID"], ascending=[False, True])

# ---------------------------
# Save concise Excel output
# ---------------------------
output_path = "polyketide_pfams_summary.xlsx"
result[["BioDolphinID", "protein_Pfam_ID", "count", "frequency_percent"]].to_excel(output_path, index=False)
print(f"Saved summary to: {output_path}")

# ---------------------------
# PIE CHART: PFAMs > 1% frequency
# ---------------------------
pfam_over1 = pfam_counts[pfam_counts["frequency_percent"] > 2].reset_index(drop=True)
if pfam_over1.empty:
    print("No PFAMs exceed 1% frequency; no pie created.")
else:
    labels = pfam_over1["protein_Pfam_ID"].tolist()
    counts = pfam_over1["count"].tolist()
    freqs = pfam_over1["frequency_percent"].tolist()

    fig, ax = plt.subplots(figsize=(7, 7))

    # Generate colors with 50% transparency
    base_colors = plt.cm.tab20(np.linspace(0, 1, len(counts)))
    base_colors[:, -1] = 0.5  # alpha = 50%
    
    wedges, texts = ax.pie(
        counts,
        labels=labels,
        startangle=90,
        textprops={'fontsize': 10},
        colors=base_colors
    )

    # Add frequency_percent text on wedges
    for i, wedge in enumerate(wedges):
        ang = (wedge.theta2 + wedge.theta1) / 2.0
        x = 0.6 * np.cos(np.deg2rad(ang))
        y = 0.6 * np.sin(np.deg2rad(ang))
        pct_text = f"{freqs[i]:.2f}%"
        txt = ax.text(x, y, pct_text, ha='center', va='center', fontsize=9, fontfamily='Arial')
        txt.set_path_effects([patheffects.Stroke(linewidth=1.5, foreground='white'), patheffects.Normal()])

    ax.set_title("PFAM distribution (>1% frequency)", fontsize=14, fontfamily='Arial')
    plt.tight_layout()
    pie_output = "polyketide_pfam_pie_chart.pdf"
    plt.savefig(pie_output, dpi=300)
    plt.close()
    print(f"Saved pie chart to: {pie_output}")
