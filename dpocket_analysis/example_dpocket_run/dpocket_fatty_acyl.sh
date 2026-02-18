#!/bin/bash

# Example SLURM script to run dpocket on a high perfomance computing cluster

# submit with: sbatch dpocket_fatty_acyl.sh

#SBATCH -J dpocket_Fatty_acyl                # Job name
#SBATCH -A gts-amcshan3                      # Account name (PACE project)
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 024:00:00                         # Wall time (adjust as needed)
#SBATCH -o logs/dpocket_%j.out               # Slurm stdout
#SBATCH -e logs/dpocket_%j.err               # Slurm stderr

# Move to the directory where sbatch was submitted
cd "$SLURM_SUBMIT_DIR"

conda activate fpocket_env

mkdir -p dpocket_outputs

FAILED_LOG="dpocket_failed.txt"
MISSING_LIGAND_PDBS="dpocket_missing_ligand_pdbs.txt"

> "$FAILED_LOG"
> "$MISSING_LIGAND_PDBS"

for input_file in split_inputs/*.txt; do
    base=$(basename "${input_file%.*}")
    output_prefix="dpocket_outputs/${base}"

    echo "Running Dpocket on $input_file..."

    # Run dpocket and capture both stdout and stderr
    dpocket_output=$(dpocket -f "$input_file" -o "$output_prefix" 2>&1)
    dpocket_exit_code=$?

    # Log general failures
    if [[ $dpocket_exit_code -ne 0 ]]; then
        echo "$input_file" >> "$FAILED_LOG"
    fi

    # Check for ligand errors and extract PDB file names
    echo "$dpocket_output" | grep "ERROR - No ligand" | while read -r line; do
        pdb_file=$(echo "$line" | grep -oP 'in \K[^.]+\.pdb')
        echo "$pdb_file" >> "$MISSING_LIGAND_PDBS"
    done
done

echo "Done!"
echo "Failed input files logged to: $FAILED_LOG"
echo "PDBs with missing ligands logged to: $MISSING_LIGAND_PDBS"
