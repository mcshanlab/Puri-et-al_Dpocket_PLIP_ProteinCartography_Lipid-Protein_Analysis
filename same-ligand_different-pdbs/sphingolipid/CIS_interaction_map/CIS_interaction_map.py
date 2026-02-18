import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Config
# -----------------------------
data_file = "CIS_lipid_atoms_total_clean.txt"
output_file = "CIS_atom_interactions_percentage.pdf"

# -----------------------------
# Load the data
# -----------------------------
df = pd.read_csv(
    data_file,
    sep=r"\s+",
    names=["AtomName", "Count", "Element"],
    dtype={"AtomName": str, "Element": str},
    na_values=["", "NA"]
)

df["Count"] = pd.to_numeric(df["Count"], errors="coerce").fillna(0)

# -----------------------------
# Sum counts per AtomName
# -----------------------------
df_grouped = (
    df.groupby(["AtomName", "Element"], as_index=False)["Count"]
    .sum()
)

# -----------------------------
# REMOVE zero-interaction atoms
# -----------------------------
df_grouped = df_grouped[df_grouped["Count"] > 0]

# If everything is zero, exit safely
if df_grouped.empty:
    print("No interacting atoms found.")
    exit()

# -----------------------------
# Calculate percentage
# -----------------------------
total = df_grouped["Count"].sum()
df_grouped["Percent"] = df_grouped["Count"] / total * 100

# Sort descending
df_grouped = df_grouped.sort_values("Percent", ascending=False)

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(10,4))

bars = plt.bar(
    df_grouped["AtomName"],
    df_grouped["Percent"],
    color="skyblue",
    edgecolor="black"
)

for bar, pct in zip(bars, df_grouped["Percent"]):
    plt.text(
        bar.get_x() + bar.get_width()/2,
        bar.get_height() + 0.5,
        f"{pct:.1f}%",
        ha="center",
        va="bottom",
        fontsize=8,
        rotation=45
    )

plt.xlabel("Atom Name")
plt.ylabel("Percentage of Interactions (%)")
plt.title("Percentage of Interactions per CIS Atom")
plt.xticks(rotation=45)
plt.ylim(0, df_grouped["Percent"].max() * 1.2)

plt.tight_layout()
plt.savefig(output_file)
plt.show()

print(f"Plot saved as {output_file}")
