#!/bin/bash

# Example SLURM script to run PLIP on a high perfomance computing cluster

# submit with: sbatch PLIP_fatty_acyl.sh

#SBATCH -J Fatty_acyl                            # Job name
#SBATCH -A gts-amcshan3                          # Charge account
#SBATCH -N 1                                     # Number of nodes
#SBATCH --ntasks-per-node=4                      # Number of tasks (cores) per node
#SBATCH --mem-per-cpu=1G                         # Memory per core
#SBATCH -t 010:00:00                             # Walltime (hh:mm:ss)
#SBATCH -q inferno                               # Queue name
#SBATCH -o Report-%j.out                         # Output file
#SBATCH --mail-type=BEGIN,END,FAIL               # Email notifications
#SBATCH --mail-user=npuri24@gatech.edu           # Replace with your email

# ----------------------------
# Job execution starts here
# ----------------------------

cd $SLURM_SUBMIT_DIR

# Activate PLIP conda environment
conda activate plip

# Define input and output folders containg PDB files (1:1 lipid pairs)
input_folder="$path/Fatty_acyl"
output_folder="$path/Fatty_acyl"

# Loop through each .pdb file in the input folder
for pdb_file in "$input_folder"/*.pdb; do
    # Full path and base filename (without .pdb)
    full_path=$(realpath "$pdb_file")
    base_name=$(basename "$pdb_file" .pdb)

    # Output folder for this specific PDB
    out_dir="$output_folder/$base_name"
    mkdir -p "$out_dir"

    # Run PLIP with XML (-x) and Text (-t) outputs
    plipcmd.py -f "$full_path" -x -t -o "$out_dir"

    echo "Processed: $base_name"
done
