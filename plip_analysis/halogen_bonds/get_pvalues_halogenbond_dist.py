import pandas as pd
import itertools
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# --------------------------------------
# Input file (donor–acceptor distances)
# --------------------------------------
input_file = "halogenbond_dist_violin_values.tsv"
output_file = "halogenbond_dist_pairwise_significance.tsv"

value_col = "halogenbond_dist"

# --------------------------------------
# Load data
# --------------------------------------
df = pd.read_csv(input_file, sep="\t")

# Ensure numeric
df[value_col] = pd.to_numeric(df[value_col], errors="coerce")
df = df.dropna(subset=[value_col])

classes = df["class"].unique()

# --------------------------------------
# Pairwise Mann–Whitney U tests
# --------------------------------------
results = []

for a, b in itertools.combinations(classes, 2):
    group_a = df[df["class"] == a][value_col]
    group_b = df[df["class"] == b][value_col]

    stat, p = mannwhitneyu(group_a, group_b, alternative="two-sided")

    results.append({
        "group_1": a,
        "group_2": b,
        "U_statistic": stat,
        "p_value": p,
        "n1": len(group_a),
        "n2": len(group_b)
    })

results_df = pd.DataFrame(results)

# --------------------------------------
# Multiple testing correction (FDR)
# --------------------------------------
reject, pvals_corrected, _, _ = multipletests(
    results_df["p_value"],
    method="fdr_bh"
)

results_df["p_adj_FDR_BH"] = pvals_corrected
results_df["significant_FDR_0.05"] = reject

# --------------------------------------
# Save results
# --------------------------------------
results_df.to_csv(output_file, sep="\t", index=False)

print(f"Saved pairwise statistics to: {output_file}")
