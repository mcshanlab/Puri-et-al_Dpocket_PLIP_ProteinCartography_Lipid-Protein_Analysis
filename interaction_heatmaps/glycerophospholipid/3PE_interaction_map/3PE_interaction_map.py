import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Config
# -----------------------------
data_file = "3PE_lipid_atoms_total_clean.txt"
output_file = "3PE_atom_interactions_percentage.pdf"

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

# Sort descending
df_grouped = df_grouped.sort_values('Percent', ascending=False).reset_index(drop=True)

# -----------------------------
# Plot
# -----------------------------
fig, ax = plt.subplots(figsize=(14,4))

# Use numeric positions to avoid categorical tick issues
x_positions = range(len(df_grouped))

bars = ax.bar(x_positions,
              df_grouped['Percent'],
              color='skyblue',
              edgecolor='black')

# -----------------------------
# Format x-axis tick labels
# -----------------------------
ax.set_xticks(x_positions)
ax.set_xticklabels(df_grouped['AtomName'],
                   rotation=45,
                   ha='right',
                   fontsize=7,
                   fontname='Arial',
                   fontweight='bold')

# -----------------------------
# Add percentage labels above bars
# -----------------------------
for bar, pct in zip(bars, df_grouped['Percent']):
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() + 0.5,
            f"{pct:.1f}%",
            ha='center',
            va='bottom',
            fontsize=8,
            rotation=45)

# -----------------------------
# Axis labels and title
# -----------------------------
ax.set_xlabel("Atom Name", fontsize=12)
ax.set_ylabel("Percentage of Interactions (%)", fontsize=12)
ax.set_title("Percentage of Interactions per 3PE Atom", fontsize=14)

ax.set_ylim(0, max(df_grouped['Percent']) * 1.2)

plt.tight_layout()
plt.savefig(output_file)
plt.show()

print(f"Plot saved as {output_file}")
