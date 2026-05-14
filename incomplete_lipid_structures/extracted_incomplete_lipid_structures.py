import os
import requests
from Bio.PDB import PDBParser

# ============================================================
# SETTINGS
# ============================================================

ROOT_DIR = "/Volumes/GigiMurin/plip"

# Ignore hydrogens?
IGNORE_HYDROGENS = True

# Save report
OUTPUT_FILE = "ligand_missing_atom_report.tsv"

# ============================================================
# CCD CACHE
# ============================================================

ccd_cache = {}

# ============================================================
# FETCH CCD ATOMS
# ============================================================

def fetch_ccd_atoms(ccd_id):
    """
    Download CCD definition from RCSB and extract atom names.
    """

    ccd_id = ccd_id.upper()

    if ccd_id in ccd_cache:
        return ccd_cache[ccd_id]

    print(f"Downloading CCD definition for: {ccd_id}")

    url = f"https://files.rcsb.org/ligands/view/{ccd_id}.cif"

    try:
        response = requests.get(url, timeout=20)
        response.raise_for_status()

    except Exception as e:
        print(f"ERROR downloading CCD {ccd_id}: {e}")
        return None

    lines = response.text.splitlines()

    expected_atoms = set()

    capture = False

    for line in lines:

        if "_chem_comp_atom.atom_id" in line:
            capture = True
            continue

        if capture:

            if line.startswith("#"):
                break

            fields = line.split()

            if len(fields) < 2:
                continue

            atom_name = fields[1].strip()

            # Ignore hydrogens
            if IGNORE_HYDROGENS and atom_name.startswith("H"):
                continue

            expected_atoms.add(atom_name)

    ccd_cache[ccd_id] = expected_atoms

    return expected_atoms

# ============================================================
# EXTRACT CCD ID DIRECTLY FROM PDB FILE
# ============================================================

def extract_ccd_from_pdb(pdb_file):
    """
    Extract ligand CCD ID directly from HETATM records.
    """

    try:

        with open(pdb_file, "r", errors="ignore") as f:

            for line in f:

                if line.startswith("HETATM"):

                    resname = line[17:20].strip()

                    # Skip common non-ligands
                    if resname in ["HOH", "WAT"]:
                        continue

                    return resname

    except Exception as e:
        print(f"ERROR reading {pdb_file}: {e}")

    return None

# ============================================================
# EXTRACT OBSERVED ATOMS FROM PDB
# ============================================================

def extract_observed_atoms(pdb_file, ccd_id):
    """
    Extract ligand atom names from HETATM records.
    """

    parser = PDBParser(QUIET=True)

    try:
        structure = parser.get_structure("x", pdb_file)

    except Exception as e:
        print(f"ERROR parsing {pdb_file}: {e}")
        return None

    observed_atoms = set()

    for model in structure:

        for chain in model:

            for residue in chain:

                if residue.resname.strip() != ccd_id:
                    continue

                for atom in residue:

                    atom_name = atom.get_name().strip()

                    # Ignore hydrogens
                    if IGNORE_HYDROGENS and atom_name.startswith("H"):
                        continue

                    observed_atoms.add(atom_name)

    return observed_atoms

# ============================================================
# MAIN ANALYSIS
# ============================================================

results = []

total_files = 0

for dirpath, dirnames, filenames in os.walk(ROOT_DIR):

    # Skip macOS hidden files
    pdb_files = [
        f for f in filenames
        if f.endswith(".pdb") and not f.startswith("._")
    ]

    if not pdb_files:
        continue

    for pdb_file in pdb_files:

        total_files += 1

        full_path = os.path.join(dirpath, pdb_file)

        print("\n===================================================")
        print(f"Analyzing file #{total_files}")
        print(f"PDB file: {full_path}")

        # ----------------------------------------------------
        # Extract CCD directly from PDB
        # ----------------------------------------------------

        ccd_id = extract_ccd_from_pdb(full_path)

        if not ccd_id:
            print("WARNING: Could not determine CCD ID")
            continue

        print(f"CCD ID: {ccd_id}")

        # ----------------------------------------------------
        # Fetch CCD atom definitions
        # ----------------------------------------------------

        expected_atoms = fetch_ccd_atoms(ccd_id)

        if expected_atoms is None:
            continue

        # ----------------------------------------------------
        # Extract observed atoms
        # ----------------------------------------------------

        observed_atoms = extract_observed_atoms(full_path, ccd_id)

        if observed_atoms is None:
            continue

        # ----------------------------------------------------
        # Compare atoms
        # ----------------------------------------------------

        missing_atoms = sorted(expected_atoms - observed_atoms)
        extra_atoms = sorted(observed_atoms - expected_atoms)

        print(f"Expected atoms: {len(expected_atoms)}")
        print(f"Observed atoms: {len(observed_atoms)}")
        print(f"Missing atoms: {len(missing_atoms)}")
        print(f"Extra atoms: {len(extra_atoms)}")

        if missing_atoms:
            print(f"Missing atoms list: {', '.join(missing_atoms)}")

        if extra_atoms:
            print(f"Extra atoms list: {', '.join(extra_atoms)}")

        # ----------------------------------------------------
        # Save results
        # ----------------------------------------------------

        results.append({
            "pdb_file": full_path,
            "ccd_id": ccd_id,
            "expected_atom_count": len(expected_atoms),
            "observed_atom_count": len(observed_atoms),
            "missing_atom_count": len(missing_atoms),
            "extra_atom_count": len(extra_atoms),
            "missing_atoms": ",".join(missing_atoms),
            "extra_atoms": ",".join(extra_atoms),
        })

# ============================================================
# WRITE REPORT
# ============================================================

with open(OUTPUT_FILE, "w") as f:

    header = [
        "pdb_file",
        "ccd_id",
        "expected_atom_count",
        "observed_atom_count",
        "missing_atom_count",
        "extra_atom_count",
        "missing_atoms",
        "extra_atoms",
    ]

    f.write("\t".join(header) + "\n")

    for r in results:

        row = [str(r[h]) for h in header]

        f.write("\t".join(row) + "\n")

# ============================================================
# DONE
# ============================================================

print("\n===================================================")
print("Analysis complete.")
print(f"Report written to: {OUTPUT_FILE}")
print(f"Structures analyzed: {len(results)}")
print("===================================================")
