import pandas as pd
import numpy as np
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

go_df = go_df[
    ["BioDolphinID", "protein_Cellular_Component_(GO)", "protein_Organism"]
].copy()

go_df["PDB_ID"] = go_df["BioDolphinID"].apply(clean_biodolphin_id)
go_df = go_df.dropna(subset=["PDB_ID"])

# --------------------------------------
# EXTRACT GO IDS
# --------------------------------------
def extract_go_ids(x):
    if pd.isna(x):
        return []
    return re.findall(r"GO:\d{7}", str(x))

go_df["GO_IDs"] = go_df["protein_Cellular_Component_(GO)"].apply(extract_go_ids)
go_df = go_df.explode("GO_IDs")
go_df = go_df.dropna(subset=["GO_IDs"])

# --------------------------------------
# TARGET GO TERMS
# --------------------------------------
targets = {
    "Saccharolipid": {"GO:0009279", "GO:0009427"},
    "Glycerolipid": {"GO:0005751", "GO:0045277"}
}

# --------------------------------------
# MAIN LOOP
# --------------------------------------
for sheet in ["Saccharolipid", "Glycerolipid"]:

    df = pd.read_excel(excel_file, sheet_name=sheet)

    pock_col = [c for c in df.columns if "pock" in c.lower()][0]
    pdb_col = "pdb" if "pdb" in df.columns else df.columns[0]

    tmp = df[[pdb_col, pock_col]].copy()
    tmp.columns = ["PDB_ID", "lig_vol"]

    tmp["PDB_ID"] = tmp["PDB_ID"].apply(clean_dpocket_id)
    tmp["lig_vol"] = pd.to_numeric(tmp["lig_vol"], errors="coerce")
    tmp = tmp.dropna()

    merged = pd.merge(tmp, go_df, on="PDB_ID", how="inner")

    print("\n" + "=" * 80)
    print(f"{sheet} — HIGH LIGAND VOLUME CONTRIBUTORS")
    print("=" * 80)

    go_set = targets[sheet]
    sub = merged[merged["GO_IDs"].isin(go_set)]

    if sub.empty:
        print("No matches found\n")
        continue

    # sort by HIGH ligand volume (this is the key fix)
    sub = sub.sort_values("lig_vol", ascending=False)

    for _, r in sub.head(50).iterrows():
        print(
            f"{r['PDB_ID']}\t"
            f"{r['GO_IDs']}\t"
            f"{r['lig_vol']:.3f}\t"
            f"{r['protein_Organism']}"
        )
