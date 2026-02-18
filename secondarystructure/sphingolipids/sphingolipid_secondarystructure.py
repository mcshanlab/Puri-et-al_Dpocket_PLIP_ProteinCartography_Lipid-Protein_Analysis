#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
from collections import defaultdict
from Bio.PDB import PDBParser, DSSP

# ============================================================
# USER SETTINGS
# ============================================================
ROOT_DIR = "/Volumes/GigiMurin/plip/sphingo_lipids"
METADATA_TSV = (
    "/Users/amcshan3/Desktop/Manuscripts/PLIP_Dpocket_Lipid_Puri_2025/used-to-make/proteincartography/final_results_cluster-mode_sphingolipids/cluster-mode-sphingolipids_aggregated_features_pca_tsne.tsv"
)
DIST_CUTOFF = 5.0  # Å

OUTPUT_PER_STRUCT = "ss_per_structure.csv"
OUTPUT_PER_LC = "ss_by_LeidenCluster.csv"

# ============================================================
# LOAD METADATA
# ============================================================
print("Loading metadata...")
meta = pd.read_csv(METADATA_TSV, sep="\t")
required_cols = {"protid", "LeidenCluster"}
missing = required_cols - set(meta.columns)
if missing:
    raise ValueError(f"Missing required metadata columns: {missing}")

protid_to_lc = dict(zip(meta["protid"], meta["LeidenCluster"]))
meta_protids = set(meta["protid"])

# ============================================================
# PDB PARSER
# ============================================================
parser = PDBParser(QUIET=True)
rows = []

# ============================================================
# WALK PDB FILES
# ============================================================
for root, _, files in os.walk(ROOT_DIR):
    for file in files:
        if not file.lower().endswith(".pdb"):
            continue

        protid_raw = os.path.splitext(file)[0]

        # Remove known prefixes/suffixes
        for prefix in ["plipfixed."]:
            if protid_raw.startswith(prefix):
                protid_raw = protid_raw[len(prefix):]
        for suffix in ["_protonated"]:
            if protid_raw.endswith(suffix):
                protid_raw = protid_raw[:-len(suffix)]

        # Relaxed match: pick the longest substring in metadata
        protid_candidates = [p for p in meta_protids if p in protid_raw]
        if protid_candidates:
            protid = max(protid_candidates, key=len)
        else:
            print("Skipping, protid not in metadata:", protid_raw)
            continue

        pdb_path = os.path.join(root, file)
        print(f"Processing: {protid}")

        try:
            structure = parser.get_structure(protid, pdb_path)
            model = structure[0]

            # ------------------------------------------------
            # DSSP (secondary structure)
            # ------------------------------------------------
            dssp = DSSP(model, pdb_path)
            residue_ss = {}
            for key in dssp.keys():
                chain_id, res_id = key
                ss = dssp[key][2]
                residue_ss[(chain_id, res_id[1])] = ss

            # ------------------------------------------------
            # Collect ligand atom coordinates (any HETATM except water)
            # ------------------------------------------------
            ligand_coords = []
            ligand_ids = set()
            for chain in model:
                for residue in chain:
                    if residue.id[0] != " " and residue.resname not in {"HOH", "WAT"}:
                        ligand_ids.add(residue.resname)
                        for atom in residue:
                            ligand_coords.append(atom.coord)

            if not ligand_coords:
                print(f"No ligands detected for {protid}, skipping...")
                continue

            ligand_coords = np.array(ligand_coords)
            print(f"Ligands found in {protid}: {ligand_ids}")

            # ------------------------------------------------
            # Find protein residues within cutoff
            # ------------------------------------------------
            ss_counts = defaultdict(set)
            for chain in model:
                for residue in chain:
                    if residue.id[0] != " ":
                        continue

                    # CORRECTED ACCESS
                    ca = residue["CA"] if "CA" in residue else None
                    if ca is None:
                        continue

                    dists = np.linalg.norm(ligand_coords - ca.coord, axis=1)
                    if np.any(dists <= DIST_CUTOFF):
                        ss = residue_ss.get((chain.id, residue.id[1]), "-")
                        ss_counts[ss].add((chain.id, residue.id[1]))

            helix = len(ss_counts.get("H", []))
            sheet = len(ss_counts.get("E", []))
            loop = sum(len(v) for k, v in ss_counts.items() if k not in ("H", "E"))
            total = helix + sheet + loop

            if total == 0:
                print(f"No residues near ligand within cutoff for {protid}, skipping...")
                continue

            rows.append({
                "protid": protid,
                "LeidenCluster": protid_to_lc[protid],
                "helix": helix,
                "sheet": sheet,
                "loop": loop,
                "total": total,
                "helix_pct": 100 * helix / total,
                "sheet_pct": 100 * sheet / total,
                "loop_pct": 100 * loop / total,
            })

        except Exception as e:
            print(f"FAILED: {protid} ({e})")

# ============================================================
# BUILD DATAFRAME AND SAVE
# ============================================================
ss_df = pd.DataFrame(rows)
if ss_df.empty:
    raise RuntimeError(
        "No structures passed all filters. "
        "Check protid matching, ligand detection, and DSSP."
    )

ss_df.to_csv(OUTPUT_PER_STRUCT, index=False)
print(f"Per-structure CSV written: {OUTPUT_PER_STRUCT}")

# ============================================================
# AGGREGATE BY LEIDEN CLUSTER
# ============================================================
lc_summary = (
    ss_df
    .groupby("LeidenCluster")
    .agg(
        n_structures=("protid", "count"),
        helix_mean=("helix_pct", "mean"),
        helix_std=("helix_pct", "std"),
        sheet_mean=("sheet_pct", "mean"),
        sheet_std=("sheet_pct", "std"),
        loop_mean=("loop_pct", "mean"),
        loop_std=("loop_pct", "std"),
    )
    .reset_index()
)

lc_summary.to_csv(OUTPUT_PER_LC, index=False)
print(f"Per-LC summary CSV written: {OUTPUT_PER_LC}")
print("\nDone.")
