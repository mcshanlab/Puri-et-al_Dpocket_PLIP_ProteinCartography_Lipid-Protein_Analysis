import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import ast

# --------------------------------------
# Config
# --------------------------------------
excel_file = "../../Source_Data_dpocket.xlsx"
bio_file = "../BioDolphin_vr1.1.csv"

target_sheets = [
    "Sterol", "Polyketide", "Prenol", "Saccharolipid",
    "Sphingolipid", "Fatty Acyl", "Glycerophospholipid", "Glycerolipid"
]

lig_col_keyword = "lig_vol"
pock_col_keyword = "pock_vol"

# --------------------------------------
# COLOR MAP
# --------------------------------------
color_map = {
    "Sterol": "#F2C9D1",
    "Polyketide": "#CBC7D6",
    "Prenol": "#D0F2F2",
    "Saccharolipid": "#F5F5F5",
    "Sphingolipid": "#F2ECD3",
    "Fatty Acyl": "#F2DDBF",
    "Glycerophospholipid": "#F5C1CE",
    "Glycerolipid": "#D4EBD1"
}

# --------------------------------------
# CLEAN IDS
# --------------------------------------
def clean_dpocket_id(x):
    if pd.isna(x):
        return np.nan
    return str(x).strip().lower().replace(".pdb", "")

def clean_biodolphin_id(x):
    if pd.isna(x):
        return np.nan
    return str(x).strip().lower()

# --------------------------------------
# Load BioDolphin
# --------------------------------------
go_df = pd.read_csv(bio_file, low_memory=False)

go_df = go_df[[
    "BioDolphinID",
    "protein_Cellular_Component_(GO)"
]].copy()

go_df["PDB_ID"] = go_df["BioDolphinID"].apply(clean_biodolphin_id)
go_df = go_df.dropna(subset=["PDB_ID"])

# --------------------------------------
# GO parsing
# --------------------------------------
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

# =========================================================
# MAIN LOOP
# =========================================================
for sheet in target_sheets:

    df = pd.read_excel(excel_file, sheet_name=sheet)

    # find columns
    lig_col = None
    pock_col = None

    for col in df.columns:
        if lig_col_keyword in col.lower():
            lig_col = col
        if pock_col_keyword in col.lower():
            pock_col = col

    if lig_col is None or pock_col is None:
        print(f"WARNING: missing columns in {sheet}")
        continue

    pdb_col = "pdb" if "pdb" in df.columns else df.columns[0]

    tmp = df[[pdb_col, lig_col, pock_col]].copy()
    tmp.columns = ["PDB_ID", "lig_vol", "pock_vol"]

    tmp["PDB_ID"] = tmp["PDB_ID"].apply(clean_dpocket_id)

    tmp["lig_vol"] = pd.to_numeric(
        tmp["lig_vol"], errors="coerce"
    )
    tmp["pock_vol"] = pd.to_numeric(
        tmp["pock_vol"], errors="coerce"
    )

    tmp = tmp.dropna(
        subset=["PDB_ID", "lig_vol", "pock_vol"]
    )

    # avoid divide by zero
    tmp = tmp[tmp["pock_vol"] > 0.000000001]

    # calculate ratio
    tmp["ratio"] = tmp["lig_vol"] / tmp["pock_vol"]

    print(f"\n================ {sheet} ================")
    print("dpocket PDBs:", tmp["PDB_ID"].nunique())
    print("BioDolphin PDBs:", go_df["PDB_ID"].nunique())

    merged = pd.merge(
        tmp,
        go_df,
        on="PDB_ID",
        how="inner"
    )

    if len(merged) == 0:
        continue

    merged["class"] = sheet

    # =====================================================
    # SAVE TSV
    # =====================================================
    out_tsv = (
        f"lipid_to_pocket_ratio_{sheet}_GO_CC_values.tsv"
    )

    merged[
        ["class", "GO_CC", "ratio"]
    ].to_csv(
        out_tsv,
        sep="\t",
        index=False
    )

    print(f"saved TSV: {out_tsv}")

    # =====================================================
    # GROUP BY GO
    # =====================================================
    groups = merged.groupby("GO_CC")["ratio"]

    top_groups = (
        groups.size()
        .sort_values(ascending=False)
        .head(10)
        .index
    )

    data = []
    labels = []
    sample_sizes = []

    for g in top_groups:
        vals = groups.get_group(g)

        data.append(vals)
        labels.append(g)
        sample_sizes.append(len(vals))

    if len(data) == 0:
        continue

    # =====================================================
    # PLOT
    # =====================================================
    plt.figure(figsize=(8, 4))

    vp = plt.violinplot(
        data,
        showmeans=False,
        showextrema=False,
        showmedians=True
    )

    for body in vp["bodies"]:
        body.set_facecolor(
            color_map.get(sheet, "#A7C7E7")
        )
        body.set_edgecolor("black")
        body.set_linewidth(1)
        body.set_alpha(0.95)

    vp["cmedians"].set_color("black")
    vp["cmedians"].set_linewidth(1.5)

    plt.xticks(
        range(1, len(labels)+1),
        labels,
        rotation=40,
        ha="right"
    )

    plt.ylabel(
        "Lipid volume / Protein pocket volume"
    )

    plt.title(sheet)

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # fixed axis across all plots
    plt.yticks([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    plt.ylim(-0.05, 3.0)

    ymax = 3.0

    for i, n in enumerate(sample_sizes, start=1):
        plt.text(
            i,
            ymax * 1.03,
            f"n = {n}",
            ha="center",
            fontsize=8
        )

    plt.subplots_adjust(
        left=0.12,
        right=0.98,
        bottom=0.35,
        top=0.88
    )

    plt.savefig(
        f"{sheet}_lipid_to_pocket_ratio_by_GO_CC.pdf",
        dpi=300
    )
    plt.close()

    print(
        f"saved plot: "
        f"{sheet}_lipid_to_pocket_ratio_by_GO_CC.pdf"
    )
