# A systematic analysis of lipid-protein interactions in the Protein Data Bank

## Introduction

This repository provides code and datasets for the systematic analysis of lipid-protein interactions across eight major lipid classes using PLIP, dpocket, and ProteinCartography. 

Key functions of this repository:

PLIP analysis: Automated identification of lipid-protein interactions at atomic resolution.

dpocket analysis: Characterization of protein binding pocket geometry and physicochemical properties.

ProteinCartography: Clustering and visualization of protein folds and interaction trends across lipid classes.

Heatmap Generation: Aggregation and visualization of lipid-protein interaction hotspots.

This repository provides all scripts, datasets, and instructions necessary to reproduce the figures and analyses presented in the publication.

## Dataset availability:

All datasets used in this study are located in separate data directories.

Source_Data_PLIP.xlsx - summary of PLIP analysis

Source_Data_dpocket.xlsx - summary of dpocket analysis

Source_Data_StatisticalSignificance_MannWhitneyU_updated.xlsx - statistical analysis

dpocket_analysis: Contains dpocket-derived pocket descriptors, including volume, polarity, and hydrophobicity. An example dpocket run is provided.

dpocket_analysis_subcellular-location Contains dpocket-derived pocket descriptors, including volume, polarity, and hydrophobicity as function of GO subcellular location.

plip_analysis: Contains PLIP-derived analysis of different types of interaction descriptors. An example PLIP run is provided.

plip_analysis_subcellular: Contains PLIP-derived analysis of different types of interaction descriptors as function of GO subcellular location.

pfam_analysis: Contains analysis of protein families (PFAMs)

secondarystructure: Analaysis of protein secondary structure elements near lipid atoms

proteincartography: Contains protein cluster assignments and fold annotations used for visualization and t-SNE plots. An example ProteinCartography run is provided.

interaction_heatmaps: Interaction "heat maps" for mapping frequency of lipid atoms contacting protein atoms

incomplete_lipid_structures: Analysis of the amount of PDB entries lacking fully modeled lipid structures.

### Prerequisites

Ensure you have the following installed:

Python ≥3.9 with Matplotlib / Seaborn / Pandas / NumPy for plotting and analysis

PLIP: https://github.com/pharmai/plip

dpocket: https://github.com/Discngine/fpocket

ProteinCartography: https://github.com/Arcadia-Science/ProteinCartography

PyMOL: https://github.com/schrodinger/pymol-open-source or https://pymol.org/

## PDB files of lipid-protein pairs:

Due to space limitations we could not provide all PDB files ( >> 100 GBs) but a list of BioDolphin IDs used is available in the Source_Data_dpocket.xlsx file and in the ProteinCartography results files (i.e., cluster-mode-sterols_aggregated_features_pca_tsne.tsv).
