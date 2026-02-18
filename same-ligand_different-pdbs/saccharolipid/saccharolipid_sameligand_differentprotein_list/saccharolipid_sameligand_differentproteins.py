#!/usr/bin/env python3

import pandas as pd
from collections import defaultdict

# Path to your TSV file
tsv_file = "/Users/amcshan3/Desktop/Manuscripts/PLIP_Dpocket_Lipid_Puri_2025/used-to-make/proteincartography/final_results_cluster-mode_saccharolipids/cluster-mode-saccharolipids_aggregated_features_pca_tsne.tsv"

# Load TSV
df = pd.read_csv(tsv_file, sep="\t", dtype=str)
df = df.dropna(subset=["protid"])

def parse_ligand(protid):
    """
    Correct ligand parsing:
    Take the first three characters of the last dash-delimited field.
    
    Examples:
    BD1jpz-A-A-1401     -> 140
    BD5uxy-A-A-X901     -> X90
    BDxxxx-A-B-LIG1402  -> LIG
    """
    last_field = protid.split("-")[-1]
    return last_field[:3] if len(last_field) >= 3 else None

def protein_signature(row):
    """
    Define a unique protein by these 4 columns:
    Protein names, LeidenCluster, Pfam, InterPro
    """
    return (
        row.get("Protein names", ""),
        row.get("LeidenCluster", ""),
        row.get("Pfam", ""),
        row.get("InterPro", "")
    )

# Build ligand -> protein_signature -> protids mapping
ligand_map = defaultdict(lambda: defaultdict(list))

for _, row in df.iterrows():
    protid = row["protid"]
    ligand = parse_ligand(protid)

    if ligand is None:
        continue

    signature = protein_signature(row)
    ligand_map[ligand][signature].append(protid)

# Print ligands shared by multiple distinct proteins
print("Ligands shared by distinct proteins:\n")

for ligand, protein_dict in sorted(ligand_map.items()):
    if len(protein_dict) > 1:
        print(f"Ligand: {ligand}")
        for i, (signature, protids) in enumerate(protein_dict.items(), start=1):
            protein_name, leiden, pfam, interpro = signature
            print(f"  Protein {i}:")
            print(f"    Protein names : {protein_name}")
            print(f"    LeidenCluster : {leiden}")
            print(f"    Pfam          : {pfam}")
            print(f"    InterPro      : {interpro}")
            print(f"    Protids:")
            for p in protids:
                print(f"      {p}")
        print()
