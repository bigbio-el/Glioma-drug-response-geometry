# Glioma drug-response geometry

MATLAB code for the single-cell analysis of drug responses across recurrent malignant glioma transcriptional states.

The analysis uses transformed glioma cells from the integrated single-cell RNA-sequencing dataset released with the consensus scHPF study. Drug responses are summarized across AC-like, OPC-like, NPC-like and MES-like states using response-direction coherence (C), change in state separation (K), and mean response magnitude (M).

## Files

`gbm_drug_geometry_analysis.m`  
Main analysis script.

`inspect_glioma_loom.m`  
Helper script for checking the structure and metadata fields of the loom file.

## Data

The analysis uses the integrated glioma loom released with:

Levitin HM, Zhao W, Bruce JN, Canoll P, Sims PA.  
**Consensus scHPF Identifies Cell Type-Specific Drug Responses in Glioma by Integrating Large-Scale scRNA-seq.**  
bioRxiv (2023).  
doi: 10.1101/2023.12.05.570193

The integrated resource contains drug-perturbed human glioma slice cultures and includes the raw count matrix together with the metadata used in the study. The source study reports 19 patients, 10 treatment conditions and 52 samples. :contentReference[oaicite:0]{index=0}

The loom file is available from the link provided by the study:

https://drive.google.com/file/d/18-KInmm43wKdBX95Gq9xbuzAQwtLjgE9/view?usp=sharing

The previously published slice-culture dataset included in this resource is GSE148842:

Zhao W et al.  
**Deconvolution of cell type-specific drug responses in human tumor tissue with single-cell RNA-seq.**  
Genome Medicine 13, 82 (2021).  
doi: 10.1186/s13073-021-00894-y

## Running the analysis

Download the integrated `.loom` file from the link above.

Place

`gbm_drug_geometry_analysis.m`

and

`inspect_glioma_loom.m`

in a folder of your choice.

In MATLAB, change the Current Folder to the folder containing the scripts and run:

```matlab
gbm_drug_geometry_analysis
