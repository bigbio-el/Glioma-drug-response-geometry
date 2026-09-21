# Glioma drug-response geometry

MATLAB code for the analysis described in:

**Single-cell geometry of drug responses across malignant glioma states**

The study reanalyzes public scRNA-seq data from drug-treated glioma slice cultures and compares responses across AC-like, OPC-like, NPC-like and MES-like malignant states.

## Files

- `Glioma_drug_geometry_analysis.m` — main analysis
- `Inspect_glioma_loom.m` — loom inspection
- `results/cross_drug_unit_summary.csv` — TissueID-level cross-drug summary
- `results/cross_drug_global_inference.csv` — cross-drug inference

## Data

Integrated loom from Levitin et al. (2023):  
doi: 10.1101/2023.12.05.570193

The dataset also includes GSE148842 from Zhao et al. (2021).

## Main analysis

The code calculates:

- response-direction coherence (C)
- change in state separation (K)
- mean response magnitude (M)

TissueID is used as the biological unit.

The final analysis includes 43 matched treatment-vehicle comparisons from 19 TissueIDs.

For the cross-drug analysis, 52 drug pairs from 11 TissueIDs were summarized once per TissueID. The within-minus-between coherence difference was positive in all 11 units.

## Run

Download the loom file, place it with the MATLAB scripts, and run:

`Glioma_drug_geometry_analysis.m`

Results are written to the `results` folder.
