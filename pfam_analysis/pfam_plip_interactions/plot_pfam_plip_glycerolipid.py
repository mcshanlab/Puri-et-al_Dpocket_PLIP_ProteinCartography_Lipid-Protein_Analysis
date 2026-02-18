#!/usr/bin/env python3
"""
Aggregate PLIP interaction counts per PFAM (frequency >2%) and plot a normalized barplot.
Y-axis is normalized per PFAM to the largest interaction count within that PFAM.

Inputs:
- Excel file with BioDolphinID, protein_Pfam_ID, count, frequency_percent
- report.txt files for each BioDolphinID containing interaction data

Outputs:
- Normalized barplot of interaction counts per PFAM (PDF with editable text)
"""

import os
import ast
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------
# Config: directories and files
# ---------------------------
excel_path = "/Users/amcshan3/Desktop/Manuscripts/PLIP_Dpocket_Lipid_Puri_2025/used-to-make/pfam/pfam_frequencies/glycerolipid_pfams_summary.xlsx"
plip_base_dir = "/Volumes/8TB_McShan_Drive/nikki/project_1/plip/Gylcerolipids"
output_plot = "pfam_plip_glycerolipid_barplot_normalized_colored.pdf"  # saved in current directory

# Ensure editable PDF text
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'

# ---------------------------
# Load Excel file
# ---------------------------
df = pd.read_excel(excel_path)
required_cols = ["BioDolphinID", "protein_Pfam_ID", "frequency_percent"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns in Excel: {missing}")

# ---------------------------
# Filter PFAMs by frequency > 2%
# ---------------------------
df_filtered = df[df["frequency_percent"] > 2]

# ---------------------------
# Define interaction types and mapping for custom colors
# ---------------------------
interaction_types = [
    'hydrophobic',
    'hydrogen_bonds',
    'salt_bridges',
    'water_bridges',
    'pi_stacking',
    'pi_cation',
    'halogen_bonds',
    'metal_bonds'
]

interaction_labels = {
    'hydrophobic': "Hydrophobic Interactions",
    'hydrogen_bonds': "Hydrogen Bonds",
    'salt_bridges': "Salt Bridges",
    'water_bridges': "Water Bridges",
    'pi_stacking': "pi-Stacking",
    'pi_cation': "pi-Cation Interactions",
    'halogen_bonds': "Halogen Bonds",
    'metal_bonds': "Metal Complexes"
}

colors = {
    "Hydrophobic Interactions": "#E88FBC",
    "Hydrogen Bonds": "#7DD3F6",
    "Salt Bridges": "#FFD478",
    "Water Bridges": "#D3D3D3",          # optional gray for water bridges
    "pi-Stacking": "#90C96A",
    "pi-Cation Interactions": "#009152",
    "Halogen Bonds": "#F79420",
    "Metal Complexes": "#C2C0C0"
}

# ---------------------------
# Function to extract bond counts from report.txt
# ---------------------------
def extract_bond_entries(filepath):
    bond_counts = {k: 0 for k in interaction_types}

    with open(filepath, 'r') as f:
        lines = f.readlines()

    current_section = None
    for line in lines:
        line_lower = line.strip().lower()
        if line_lower.startswith('**'):
            if 'hydrophobic' in line_lower:
                current_section = 'hydrophobic'
            elif 'hydrogen' in line_lower:
                current_section = 'hydrogen_bonds'
            elif 'salt' in line_lower:
                current_section = 'salt_bridges'
            elif 'water' in line_lower:
                current_section = 'water_bridges'
            elif 'stacking' in line_lower:
                current_section = 'pi_stacking'
            elif 'cation' in line_lower:
                current_section = 'pi_cation'
            elif 'halogen' in line_lower:
                current_section = 'halogen_bonds'
            elif 'metal' in line_lower:
                current_section = 'metal_bonds'
            else:
                current_section = None
            continue
        if current_section and line.startswith('|') and 'RESNR' not in line:
            bond_counts[current_section] += 1

    return bond_counts

# ---------------------------
# Aggregate counts per PFAM
# ---------------------------
pfam_list = df_filtered['protein_Pfam_ID'].dropna().unique().tolist()
pfam_agg = []

for pfam_entry in pfam_list:
    # Parse pfam_entry if stored as a list string
    if isinstance(pfam_entry, str) and pfam_entry.startswith('['):
        try:
            pfams = ast.literal_eval(pfam_entry)
        except Exception:
            pfams = [pfam_entry]
    else:
        pfams = [pfam_entry]

    # Initialize counts
    total_counts = {k: 0 for k in interaction_types}

    # Find all BioDolphinIDs corresponding to this PFAM entry
    pfam_mask = df_filtered['protein_Pfam_ID'] == pfam_entry
    bio_ids = df_filtered.loc[pfam_mask, 'BioDolphinID'].dropna().tolist()

    for bd_id in bio_ids:
        report_path = os.path.join(plip_base_dir, bd_id, 'report.txt')
        if os.path.exists(report_path):
            bond_counts = extract_bond_entries(report_path)
            for k in interaction_types:
                total_counts[k] += bond_counts[k]
        else:
            print(f"⚠️ report.txt missing for {bd_id}")

    for pfam in pfams:
        pfam_agg.append({'PFAM': pfam, **total_counts})

# ---------------------------
# Create DataFrame for plotting
# ---------------------------
plot_df = pd.DataFrame(pfam_agg)

# ---------------------------
# Normalize counts per PFAM
# ---------------------------
for pfam in plot_df['PFAM'].unique():
    mask = plot_df['PFAM'] == pfam
    max_val = plot_df.loc[mask, interaction_types].values.max()
    if max_val > 0:
        plot_df.loc[mask, interaction_types] = plot_df.loc[mask, interaction_types] / max_val

# Melt for plotting and map labels
plot_df_melt = plot_df.melt(id_vars='PFAM', var_name='Interaction', value_name='Normalized_Count')
plot_df_melt['Interaction_Label'] = plot_df_melt['Interaction'].map(interaction_labels)

# ---------------------------
# Plot barplot
# ---------------------------
plt.figure(figsize=(12, 6))
sns.barplot(
    data=plot_df_melt,
    x='PFAM',
    y='Normalized_Count',
    hue='Interaction_Label',
    palette=colors
)
plt.xticks(rotation=45, ha='right')
plt.ylabel('Normalized Interaction Count')
plt.title('Normalized Interaction Counts per PFAM (frequency >2%)')
plt.legend(title='Interaction Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(output_plot, dpi=300)
plt.close()

print(f"Normalized barplot with custom colors and external legend saved to: {output_plot}")
