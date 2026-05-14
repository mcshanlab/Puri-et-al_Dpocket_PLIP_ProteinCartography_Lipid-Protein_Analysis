import pandas as pd
import numpy as np
import ast
import re

# --------------------------------------
# FILES
# --------------------------------------
excel_file = "../../Source_Data_dpocket.xlsx"
bio_file = "../BioDolphin_vr1.1.csv"

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
# LOAD BIO DOLPHIN
# --------------------------------------
go_df = pd.read_csv(bio_file, low_memory=False)

go_df = go_df[[
    "BioDolphinID",
    "protein_Cellular_Component_(GO)",
    "protein_Organism"
]].copy()

go_df["PDB_ID"] = go_df["BioDolphinID"].apply(clean_biodolphin_id)
go_df = go_df.dropna(subset=["PDB_ID"])

# --------------------------------------
# EXTRACT GO IDs (IMPORTANT FIX)
# --------------------------------------
def extract_go_ids(x):
    if pd.isna(x):
        return []
    if isinstance(x, list):
        x = str(x)

    return re.findall(r"GO:\d{7}", str(x))

go_df["GO_IDs"] = go_df["protein_Cellular_Component_(GO)"].apply(extract_go_ids)
go_df = go_df.explode("GO_IDs").dropna(subset=["GO_IDs"])

# --------------------------------------
# TARGET GO SETS (WHAT YOU ASKED)
# --------------------------------------
targets = {
    "Saccharolipid_LOW": {"GO:0009279", "GO:0009427"},
    "Glycerolipid_HIGH": {"GO:0005751", "GO:0045277"}
}

# --------------------------------------
# PROCESS EACH SHEET
# --------------------------------------
for sheet in ["Saccharolipid", "Glycerolipid"]:

    df = pd.read_excel(excel_file, sheet_name=sheet)

    pock_col = [c for c in df.columns if "pock" in c.lower()][0]
    pdb_col = "pdb" if "pdb" in df.columns else df.columns[0]

    tmp = df[[pdb_col, pock_col]].copy()
    tmp.columns = ["PDB_ID", "pock_vol"]

    tmp["PDB_ID"] = tmp["PDB_ID"].apply(clean_dpocket_id)
    tmp["pock_vol"] = pd.to_numeric(tmp["pock_vol"], errors="coerce")
    tmp = tmp.dropna()

    merged = pd.merge(tmp, go_df, on="PDB_ID", how="inner")

    print("\n" + "="*80)
    print(f"{sheet} analysis")
    print("="*80)

    # --------------------------------------
    # SACCHAROLIPID (LOW VOLUME)
    # --------------------------------------
    if sheet == "Saccharolipid":

        go_set = targets["Saccharolipid_LOW"]
        sub = merged[merged["GO_IDs"].isin(go_set)]

        print("\n📉 SACCHAROLIPID LOW VOLUME CONTRIBUTORS")
        print("GO terms:", go_set)
        print("-"*80)

        if sub.empty:
            print("No matches found\n")
        else:
            for _, r in sub.sort_values("pock_vol").head(50).iterrows():
                print(f"{r['PDB_ID']}\t{r['GO_IDs']}\t{r['pock_vol']:.2f}\t{r['protein_Organism']}")

    # --------------------------------------
    # GLYCEROLIPID (HIGH VOLUME)
    # --------------------------------------
    if sheet == "Glycerolipid":

        go_set = targets["Glycerolipid_HIGH"]
        sub = merged[merged["GO_IDs"].isin(go_set)]

        print("\n📈 GLYCEROLIPID HIGH VOLUME CONTRIBUTORS")
        print("GO terms:", go_set)
        print("-"*80)

        if sub.empty:
            print("No matches found\n")
        else:
            for _, r in sub.sort_values("pock_vol", ascending=False).head(50).iterrows():
                print(f"{r['PDB_ID']}\t{r['GO_IDs']}\t{r['pock_vol']:.2f}\t{r['protein_Organism']}")
