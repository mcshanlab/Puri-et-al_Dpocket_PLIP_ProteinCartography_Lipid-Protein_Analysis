
cd ~/Desktop/ProteinCartography-main/demo/cluster-mode_sterols/input

(1) Step 1 - Copy BioDolphin PDBs

Edit copy_pdbs.py

python3 copy_pdbs.py

(2) python3 rename.py

(4) python3 generate_tsv.py

(5) python3 update_tsv.py

Delete the template uniprot_features.tsv

Rename the updated tsv file to uniprot_features.tsv

cd ~/Desktop/ProteinCartography-main/

snakemake --configfile /PATH/cluster-mode_sterols/config.yml --cores 8