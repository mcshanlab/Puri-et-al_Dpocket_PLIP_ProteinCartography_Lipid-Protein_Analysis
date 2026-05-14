import pandas as pd
import itertools
import glob
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# --------------------------------------
# INPUT: all per-class ratio TSV files
# --------------------------------------
input_files = glob.glob("lipid_to_pocket_ratio_*_GO_CC_values.tsv")

output_file = "lipid_to_pocket_ratio_GO_within_lipid_significance.tsv"

all_results = []

# --------------------------------------
# LOOP OVER EACH LIPID CLASS FILE
# --------------------------------------
for file in input_files:

    df = pd.read_csv(file, sep="\t")

    # use ratio column
    df["ratio"] = pd.to_numeric(df["ratio"], errors="coerce")
    df = df.dropna(subset=["ratio", "GO_CC", "class"])

    lipid_class = df["class"].iloc[0]

    print(f"\n================ {lipid_class} ================")
    print("rows:", len(df))

    go_counts = df["GO_CC"].value_counts()

    # remove very small groups
    go_groups = go_counts[go_counts >= 3].index.tolist()

    if len(go_groups) < 2:
        print("Skipping (not enough GO groups)")
        continue

    class_records = []
    class_pvals = []

    # --------------------------------------
    # pairwise Mann–Whitney tests
    # --------------------------------------
    for a, b in itertools.combinations(go_groups, 2):

        group_a = df[df["GO_CC"] == a]["ratio"]
        group_b = df[df["GO_CC"] == b]["ratio"]

        if len(group_a) < 2 or len(group_b) < 2:
            continue

        stat, p = mannwhitneyu(
            group_a,
            group_b,
            alternative="two-sided"
        )

        class_records.append({
            "lipid_class": lipid_class,
            "group_1": a,
            "group_2": b,
            "U_statistic": stat,
            "p_value": p,
            "n1": len(group_a),
            "n2": len(group_b)
        })

        class_pvals.append(p)

    if len(class_records) == 0:
        continue

    # --------------------------------------
    # FDR correction within each lipid class
    # --------------------------------------
    reject, p_adj, _, _ = multipletests(
        class_pvals,
        method="fdr_bh"
    )

    for rec, padj, sig in zip(class_records, p_adj, reject):
        rec["p_adj_FDR_BH"] = padj
        rec["significant_FDR_0.05"] = sig
        all_results.append(rec)

# --------------------------------------
# SAVE
# --------------------------------------
results_df = pd.DataFrame(all_results)
results_df.to_csv(output_file, sep="\t", index=False)

print(f"\nSaved: {output_file}")
