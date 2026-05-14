import pandas as pd
import itertools
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# --------------------------------------
# Input files mapped to lipid classes
# --------------------------------------
files = {
    "sterol": "./sterols/sterol_plip_interaction_counts.tsv",
    "sphingolipid": "./sphingolipids/sphingolipids_plip_interaction_counts.tsv",
    "saccharolipid": "./saccharolipids/saccharolipids_plip_interaction_counts.tsv",
    "prenol": "./prenols/prenols_plip_interaction_counts.tsv",
    "polyketide": "./polyketides/polyketides_plip_interaction_counts.tsv",
    "glycerophospholipid": "./glycerophospholipids/glycerophospholipids_plip_interaction_counts.tsv",
    "glycerolipid": "./glycerolipids/glycerolipids_plip_interaction_counts.tsv",
    "fattyacyl": "./fattyacyls/fattyacyls_plip_interaction_counts.tsv",
}

# --------------------------------------
# Load + merge
# --------------------------------------
dfs = []
for lipid_class, file in files.items():
    df_tmp = pd.read_csv(file, sep="\t")
    df_tmp["class"] = lipid_class
    dfs.append(df_tmp)

df = pd.concat(dfs, ignore_index=True)

# --------------------------------------
# Interaction types
# --------------------------------------
interaction_cols = [
    "hydrophobic",
    "hbond",
    "saltbridge",
    "pistacking",
    "pication",
    "halogen",
    "metal"
]

classes = df["class"].unique()

# --------------------------------------
# Run analysis per interaction type
# --------------------------------------
all_results = []

for interaction in interaction_cols:

    # ensure numeric safety
    df[interaction] = pd.to_numeric(df[interaction], errors="coerce")
    df_i = df.dropna(subset=[interaction])

    results = []

    for a, b in itertools.combinations(classes, 2):

        group_a = df_i[df_i["class"] == a][interaction]
        group_b = df_i[df_i["class"] == b][interaction]

        if len(group_a) == 0 or len(group_b) == 0:
            continue

        stat, p = mannwhitneyu(group_a, group_b, alternative="two-sided")

        results.append({
            "interaction_type": interaction,
            "group_1": a,
            "group_2": b,
            "U_statistic": stat,
            "p_value": p,
            "n1": len(group_a),
            "n2": len(group_b)
        })

    results_df = pd.DataFrame(results)

    # --------------------------------------
    # FDR correction (within interaction type)
    # --------------------------------------
    if len(results_df) > 0:
        reject, pvals_corr, _, _ = multipletests(
            results_df["p_value"],
            method="fdr_bh"
        )

        results_df["p_adj_FDR_BH"] = pvals_corr
        results_df["significant_FDR_0.05"] = reject

    all_results.append(results_df)

# --------------------------------------
# Combine + save
# --------------------------------------
final_df = pd.concat(all_results, ignore_index=True)

final_df.to_csv(
    "lipid_class_interaction_type_pairwise_stats.tsv",
    sep="\t",
    index=False
)

print("Saved: lipid_class_interaction_type_pairwise_stats.tsv")
