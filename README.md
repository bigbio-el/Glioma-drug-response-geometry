# Single-cell geometry of drug responses across malignant glioma states

Code and result tables for the manuscript:

**Single-cell geometry of drug responses across malignant glioma states**  
Evgeniia Lavrenteva and Kyung Do Kim

This repository contains MATLAB code for a reanalysis of drug-treated patient-derived glioma slice-culture scRNA-seq data. The analysis asks whether AC-like, OPC-like, NPC-like and MES-like malignant programs respond in similar transcriptional directions under the same perturbation.

The main quantities are:

- **C** — response-direction coherence, defined from pairwise cosine similarity among state-conditioned treatment vectors.
- **K** — log2 treated/vehicle change in mean pairwise state-centroid separation.
- **M** — mean magnitude of the state-conditioned treatment vectors.

The analysis uses **TissueID as the biological unit**.

## Dataset

The integrated resource was released with:

Levitin HM, Zhao W, Bruce JN, Canoll P, Sims PA.  
**Consensus scHPF identifies cell type-specific drug responses in glioma by integrating large-scale scRNA-seq.**  
bioRxiv (2023). doi: 10.1101/2023.12.05.570193

It contains 391,444 cells. The present analysis retains 252,106 cells annotated by the source study as transformed glioma cells and evaluates 43 matched treatment-vehicle comparisons from 19 TissueID-defined biological units.

The integrated loom file is available from the source study:
https://drive.google.com/file/d/18-KInmm43wKdBX95Gq9xbuzAQwtLjgE9/view?usp=sharing

The previously published slice-culture dataset included in the resource is GSE148842:

Zhao W et al.  
**Deconvolution of cell type-specific drug responses in human tumor tissue with single-cell RNA-seq.**  
Genome Medicine 13, 82 (2021). doi: 10.1186/s13073-021-00894-y

## Repository contents

- `Glioma_drug_geometry_analysis.m` — main analysis script currently in the repository.
- `Inspect_glioma_loom.m` — helper script for inspecting loom structure and metadata.
- `results/cross_drug_unit_summary.csv` — TissueID-level cross-drug specificity summaries.
- `results/cross_drug_global_inference.csv` — global TissueID-level inference for the cross-drug analysis.

## Main reported results

Across the 43 treatment-vehicle comparisons, response-direction coherence was positive throughout and remained concordant under multiple representation and state-definition sensitivity analyses.

For the cross-drug specificity analysis, 52 drug pairs were available from 11 biological units. Because these pairs are nested within TissueID, inference was performed after reducing the contrasts to one median difference per TissueID.

For within-drug coherence minus between-drug coherence across different states:

- 11 biological units
- 11/11 positive differences
- median difference = 0.291
- 95% bootstrap CI = 0.201 to 0.495
- one-sided exact sign-test P = 0.000488

For the complementary same-state comparison:

- 9/11 positive differences
- median difference = 0.196
- 95% bootstrap CI = 0.057 to 0.430
- one-sided exact sign-test P = 0.0327

These analyses describe transcriptional response geometry and are not drug-efficacy scores.

## Running the analysis

1. Download the integrated `.loom` file from the source-study link above.
2. Place the MATLAB scripts in the same working folder.
3. Run the loom-inspection helper if needed.
4. Run the main analysis script and select the loom file when prompted.
5. Output tables are written to a `results` directory.

The analysis uses fixed inclusion rules and treats cells as observations within biological units rather than as independent biological replicates.

## Reproducibility

The manuscript reports sensitivity analyses based on high-confidence state filtering, an independently selected 500-gene space, published Neftel modules, non-overlapping gene cross-fitting, within-unit null calibration, and TissueID-level cross-drug specificity analysis.

The repository is intended to accompany the manuscript and to provide the code and machine-readable summaries needed to reproduce the reported computational analyses.

## Citation

Please cite the manuscript associated with this repository once bibliographic details are available.
