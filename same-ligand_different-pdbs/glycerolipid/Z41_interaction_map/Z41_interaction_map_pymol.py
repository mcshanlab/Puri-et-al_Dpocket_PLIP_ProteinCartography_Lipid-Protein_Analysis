from pymol import cmd

# Load PDB
cmd.load("Z41.pdb", "Z41")
cmd.show("sticks", "Z41")

# Read atom counts file
atom_counts = {}
total = 0
with open("Z41_lipid_atoms_total_clean.txt") as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 2:
            name = parts[0]
            count = int(parts[1])
            atom_counts[name] = count
            total += count

# Compute percentages and assign colors
max_count = max(atom_counts.values())
for atom_name, count in atom_counts.items():
    # simple red-to-blue scale
    ratio = count / max_count
    r = ratio
    g = 0
    b = 1 - ratio
    cmd.set_color(f"color_{atom_name}", [r, g, b])
    cmd.color(f"color_{atom_name}", f"Z41 and name {atom_name}")

cmd.zoom("Z41")
