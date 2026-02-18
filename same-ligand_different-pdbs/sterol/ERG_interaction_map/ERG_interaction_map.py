import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Config
# -----------------------------
data_file = "ERG_lipid_atoms_total_clean.txt"
output_file = "ERG_atom_interactions_percentage.pdf"

# -----------------------------
# Load the data
# -----------------------------
# Columns: AtomName, Count, Element
df = pd.read_csv(data_file, sep="\s+", names=["AtomName", "Count", "Element"])

# -----------------------------
# Sum counts per AtomName
# -----------------------------
df_grouped = df.groupby(['AtomName', 'Element'], as_index=False)['Count'].sum()

# -----------------------------
# Calculate percentage
# -----------------------------
df_grouped['Percent'] = df_grouped['Count'] / df_grouped['Count'].sum() * 100

# Sort descending so hottest atoms appear first
df_grouped = df_grouped.sort_values('Percent', ascending=False)

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(10,4))
bars = plt.bar(df_grouped['AtomName'], df_grouped['Percent'], color='skyblue', edgecolor='black')

# Add percentage labels above bars
for bar, pct in zip(bars, df_grouped['Percent']):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
             f"{pct:.1f}%", ha='center', va='bottom', fontsize=8)

plt.xlabel("Atom Name")
plt.ylabel("Percentage of Interactions (%)")
plt.title("Percentage of Interactions per ERG Atom")
plt.xticks(rotation=45)
plt.ylim(0, max(df_grouped['Percent'])*1.2)

plt.tight_layout()
plt.savefig(output_file)
plt.show()

print(f"Plot saved as {output_file}")
