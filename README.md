# Glioma drug-response geometry

MATLAB code used for the single-cell analysis of drug responses across recurrent malignant glioma transcriptional states.

The analysis uses the integrated glioma loom released with the consensus scHPF study by Levitin et al. and focuses on transformed glioma cells. Drug responses are summarized across AC-like, OPC-like, NPC-like and MES-like states using response-direction coherence (C), change in state separation (K), and mean response magnitude (M).

## Files

`gbm_drug_geometry_analysis.m`  
Main analysis script.

`inspect_glioma_loom.m`  
Small helper script for checking the structure and metadata fields of the loom file.

## Data

The integrated loom used in the analysis is available from the consensus scHPF study:

https://drive.google.com/file/d/18-KInmm43wKdBX95Gq9xbuzAQwtLjgE9/view?usp=sharing

The loom contains the count matrix and metadata for the integrated drug-perturbation dataset.

## Running the analysis

Download the loom file and place the MATLAB scripts in a folder of your choice.

In MATLAB, set the Current Folder to the folder containing the scripts and run:

```matlab
gbm_drug_geometry_analysis
