import pandas as pd
import plotly.graph_objects as go

# ============================================================
# INPUT FILE
# ============================================================

INPUT_FILE = "ligand_missing_atom_report.tsv"

# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv(INPUT_FILE, sep="\t")

df["expected_atom_count"] = pd.to_numeric(df["expected_atom_count"])
df["observed_atom_count"] = pd.to_numeric(df["observed_atom_count"])

# ============================================================
# DEFINE COMPLETENESS (BINARY)
# ============================================================

df["complete"] = (
    df["expected_atom_count"] == df["observed_atom_count"]
)

# ============================================================
# ATOM-LEVEL MISSINGNESS
# ============================================================

df["missing_atoms"] = (
    df["expected_atom_count"] - df["observed_atom_count"]
)

df["fraction_missing_atoms"] = (
    df["missing_atoms"] / df["expected_atom_count"]
)

df["percent_atoms_missing"] = (
    100 * df["fraction_missing_atoms"]
)

# ============================================================
# DATASET-LEVEL SUMMARY
# ============================================================

total_structures = len(df)
total_ccds = df["ccd_id"].nunique()

complete_structures = df["complete"].sum()
incomplete_structures = total_structures - complete_structures

percent_complete_structures = (
    100 * complete_structures / total_structures
)

# ============================================================
# INCOMPLETE STRUCTURE SUMMARY
# ============================================================

incomplete_df = df[~df["complete"]]

avg_missing_atoms = (
    incomplete_df["percent_atoms_missing"].mean()
)

median_missing_atoms = (
    incomplete_df["percent_atoms_missing"].median()
)

max_missing_atoms = (
    incomplete_df["percent_atoms_missing"].max()
)

min_missing_atoms = (
    incomplete_df["percent_atoms_missing"].min()
)

# ============================================================
# CCD-LEVEL STATS
# ============================================================

ccd_summary = df.groupby("ccd_id")["complete"].mean()

ccd_counts = df.groupby("ccd_id").size()

# CCDs where ALL structures are complete
ccd_100 = (ccd_summary == 1.0).sum()

# CCDs where AT LEAST ONE structure is incomplete
ccd_partial = (ccd_summary < 1.0).sum()

# CCDs where NO structures are complete
ccd_none_complete = (ccd_summary == 0.0).sum()

# ============================================================
# OPTIONAL LOW-N FLAGS
# ============================================================

ccd_low_n = ccd_counts < 1
ccd_high_n = ccd_counts >= 1

# ============================================================
# BEST / WORST CCDs
# ============================================================

ccd_summary_filtered = ccd_summary[ccd_high_n]

best_ccds = (
    ccd_summary_filtered
    .sort_values(ascending=False)
    .head(5)
)

worst_ccds = (
    ccd_summary_filtered
    .sort_values(ascending=True)
    .head(5)
)

# ============================================================
# PRINT SUMMARY
# ============================================================

print("\n==================== DATASET SUMMARY ====================\n")

print(f"Total structures analyzed: {total_structures}")
print(f"Unique CCD IDs: {total_ccds}\n")

print(
    f"Fully complete structures: "
    f"{complete_structures} "
    f"({percent_complete_structures:.2f}%)"
)

print(
    f"Incomplete structures: "
    f"{incomplete_structures} "
    f"({100 - percent_complete_structures:.2f}%)\n"
)

print("==================== ATOM-LEVEL LOSS ====================\n")

print(
    f"Avg % atoms missing "
    f"(incomplete structures): "
    f"{avg_missing_atoms:.2f}%"
)

print(
    f"Median % atoms missing "
    f"(incomplete structures): "
    f"{median_missing_atoms:.2f}%"
)

print(f"Max % atoms missing: {max_missing_atoms:.2f}%")
print(f"Min % atoms missing: {min_missing_atoms:.2f}%\n")

print(f"CCD IDs with 100% completeness: {ccd_100}")
print(f"CCD IDs with <100% completeness: {ccd_partial}")
print(f"CCD IDs with no fully complete structures: {ccd_none_complete}\n")

print(f"CCD IDs with n < 1: {ccd_low_n.sum()}")
print(f"CCD IDs with n ≥ 1: {ccd_high_n.sum()}\n")

print("==================== TOP CCDs ====================\n")

for ccd, val in best_ccds.items():
    print(
        f"{ccd}: "
        f"{val*100:.1f}% complete "
        f"(n={ccd_counts[ccd]})"
    )

print("\n==================== WORST CCDs ====================\n")

for ccd, val in worst_ccds.items():
    print(
        f"{ccd}: "
        f"{val*100:.1f}% complete "
        f"(n={ccd_counts[ccd]})"
    )

print("\n=========================================================\n")

# ============================================================
# CCD SUMMARY FOR PLOTTING
# ============================================================

summary = (
    df.groupby("ccd_id")
    .agg(
        percent_complete=("complete", "mean"),
        n_structures=("complete", "count")
    )
    .reset_index()
)

summary["percent_complete"] *= 100

summary = summary[
    summary["n_structures"] >= 1
]

summary = summary.sort_values(
    "percent_complete",
    ascending=False
)

# ============================================================
# INTERACTIVE PLOT
# ============================================================

fig = go.Figure()

fig.add_trace(
    go.Bar(
        x=summary["ccd_id"],
        y=summary["percent_complete"],
        text=summary["percent_complete"].round(1),
        customdata=summary["n_structures"],
        hovertemplate=(
            "CCD: %{x}<br>"
            "Completeness: %{y:.1f}%<br>"
            "N structures: %{customdata}<extra></extra>"
        ),
    )
)

fig.update_layout(
    title="Ligand Structural Completeness by CCD ID",
    xaxis_title="CCD ID",
    yaxis_title="Percent complete structures (%)",
    xaxis_tickangle=-90,
    height=600,
)

fig.show()
