## Scripts to run ProteinCartography

### Step 0: Download and install ProteinCartography as described elsewhere: https://github.com/Arcadia-Science/ProteinCartography

### Step 1: Obtain PDB files of lipid-protein complexes separated into separate folders by lipid class (i.e., put all sterol-protein complexes into a single folder). PDB files are 1:1 lipid-protein pairs. 

### Step 2: Put the following scripts from fixed_cartography_scripts folder into the main/ProteinCartography directory. These scripts contain edits needed to plot larger datasets and color nodes by PFAM values.

- assess_pdbs.py
- plot_interactive.py
- semantic_analysis.py

#### Step 3: Example folder to run is provided: example_cluster-mode_sterols. Modify as needed. Move this folder into ~/ProteinCartography-main/demo/

```

cd ~/ProteinCartography-main/demo/cluster-mode_sterols/input

# Step 3a - Copy BioDolphin PDBs

Edit copy_pdbs.py with the path to your sterol-protein PDB files

python3 copy_pdbs.py

# Step 3b - Rename the PDBs (if required)

python3 rename.py

# Step 3b - Generate a blank TSV file. TSV is file required to run ProteinCartography.

python3 generate_tsv.py

# Step 3c - Update the TSV file [requires BioDolphin_vr1.1.csv that can be downloaded on the BioDolphin website]

python3 update_tsv.py

# Step 3d

Delete the original template uniprot_features.tsv

Rename the updated tsv file to uniprot_features.tsv

An example uniprot_features.tsv is provided

# Step 3e - Run ProteinCartography

cd ~/ProteinCartography-main/

snakemake --configfile /PATH/cluster-mode_sterols/config.yml --cores 8

```

### Step 4: Expected results files are also provided. 

See final_results_cluster-mode_sterols and others.

The most important files for analysis are cluster-mode-sterols_aggregated_features_pca_tsne.html and cluster-mode-sterols_aggregated_features_pca_tsne.tsv. 

The other generated files are also very useful!
