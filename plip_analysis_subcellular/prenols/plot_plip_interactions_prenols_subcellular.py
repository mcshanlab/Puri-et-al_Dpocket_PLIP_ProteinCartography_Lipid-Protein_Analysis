import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import ast

# -----------------------------
# Load interaction data
# -----------------------------
df = pd.read_csv("prenols_plip_interaction_counts.tsv", sep="\t")

# -----------------------------
# Extract BioDolphin ID from file path
# -----------------------------
df["PDB_ID"] = df["file"].str.extract(r"(BD[^/]+)/report\.txt")[0]
df["PDB_ID"] = df["PDB_ID"].str.lower()

# -----------------------------
# Define interaction order
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
# Load GO annotation
# -----------------------------
bio_file = "../BioDolphin_vr1.1.csv"

go_df = pd.read_csv(
    bio_file,
    low_memory=False
)[[
    "BioDolphinID",
    "protein_Cellular_Component_(GO)"
]]

go_df["PDB_ID"] = go_df["BioDolphinID"].str.lower()

# -----------------------------
# Parse GO terms
# -----------------------------
def parse_go(x):
    if pd.isna(x):
        return []
    if isinstance(x, list):
        return x
    if isinstance(x, str):
        x = x.strip()

        if x.startswith("["):
            try:
                return ast.literal_eval(x)
            except:
                return []

        if ";" in x:
            return [i.strip() for i in x.split(";") if i.strip()]

        return [x]

    return []

go_df["protein_Cellular_Component_(GO)"] = (
    go_df["protein_Cellular_Component_(GO)"]
    .apply(parse_go)
)

go_df = go_df.explode("protein_Cellular_Component_(GO)")
go_df = go_df.rename(
    columns={"protein_Cellular_Component_(GO)": "GO_CC"}
)
go_df = go_df.dropna(subset=["GO_CC"])
go_df = go_df[go_df["GO_CC"] != ""]
go_df["GO_CC"] = go_df["GO_CC"].str.strip()

# -----------------------------
# Merge interaction + GO
# -----------------------------
merged = pd.merge(
    df,
    go_df,
    on="PDB_ID",
    how="inner"
)

# -----------------------------
# Identify high-contributing files
# (> mean + 2*SD within each GO)
# -----------------------------
out_lines = []

for go in merged["GO_CC"].unique():

    sub = merged[merged["GO_CC"] == go]

    out_lines.append("=" * 80)
    out_lines.append(f"GO Cellular Component: {go}")
    out_lines.append("=" * 80)

    for interaction in order:

        vals = sub[interaction]

        mu = vals.mean()
        sigma = vals.std()

        threshold = mu + 2 * sigma

        hits = sub[sub[interaction] > threshold][
            ["file", interaction]
        ].sort_values(
            by=interaction,
            ascending=False
        )

        out_lines.append("")
        out_lines.append(
            f"{interaction}: "
            f"mean={mu:.3f}, "
            f"sd={sigma:.3f}, "
            f"threshold={threshold:.3f}"
        )

        if len(hits) == 0:
            out_lines.append("  none")
        else:
            for _, row in hits.iterrows():
                out_lines.append(
                    f"  {row['file']}  "
                    f"(count={row[interaction]})"
                )

    out_lines.append("\n")

with open(
    "prenols_high_interaction_contributors_by_subcellular_location.txt",
    "w"
) as f:
    f.write("\n".join(out_lines))

print(
    "\nSaved: "
    "prenols_high_interaction_contributors_by_subcellular_location.txt"
)

# -----------------------------
# Normalize per structure
# -----------------------------
merged["total"] = merged[order].sum(axis=1)
merged = merged[merged["total"] > 0]

for col in order:
    merged[col] = merged[col] / merged["total"]

# -----------------------------
# Match OTHER SCRIPT:
# choose top 10 GO by ROW COUNT
# -----------------------------
top_go = (
    merged.groupby("GO_CC")
    .size()
    .sort_values(ascending=False)
    .head(10)
    .index
)

merged = merged[merged["GO_CC"].isin(top_go)]

# -----------------------------
# Count unique structures per GO
# (used for n in legend)
# -----------------------------
go_counts = (
    merged.groupby("GO_CC")["PDB_ID"]
    .nunique()
    .sort_values(ascending=False)
)

print("\nEntries analyzed per cellular location (unique PDBs):")
print(go_counts)

print("\nRows per cellular location (used for top 10 filter):")
print(
    merged.groupby("GO_CC")
    .size()
    .sort_values(ascending=False)
)

# -----------------------------
# Mean interaction profile
# -----------------------------
groups = (
    merged.groupby("GO_CC")[order]
    .mean()
    .loc[top_go]
)

# -----------------------------
# Normalize each GO row to max=100
# -----------------------------
groups_norm = (
    groups.div(groups.max(axis=1), axis=0) * 100
)

# -----------------------------
# Plot
# -----------------------------
fig, ax = plt.subplots(figsize=(12, 7))

x = np.arange(len(order))
width = 0.8 / len(groups_norm)

for i, go in enumerate(groups_norm.index):
    ax.bar(
        x + i * width,
        groups_norm.loc[go],
        width=width,
        label=f"{go} (n={go_counts[go]})"
    )

ax.set_xticks(
    x + width * (len(groups_norm.index)-1)/2
)

ax.set_xticklabels(
    labels,
    rotation=30,
    ha="right"
)

ax.set_ylabel("Normalized interaction profile (max = 100)")
ax.set_title(
    "Prenols–Protein Interaction Profiles by Subcellular Location"
)

# -----------------------------
# Legend
# -----------------------------
ax.legend(
    title="GO Cellular Component",
    loc="upper center",
    bbox_to_anchor=(0.5, -0.22),
    ncol=3,
    frameon=False
)

plt.tight_layout()
plt.subplots_adjust(bottom=0.30)

# -----------------------------
# Save
# -----------------------------
plt.savefig(
    "prenols_interaction_profiles_by_subcellular_location.pdf",
    format="pdf",
    bbox_inches="tight"
)

plt.show()
