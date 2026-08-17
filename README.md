# Spatial-Project
# Spatial Transcriptomics Architecture & Multi-Engine Deconvolution Pipeline

[![Seurat v5](https://img.shields.io/badge/Seurat-v5.0-blue)](https://satijalab.org/seurat/)
[![spacexr / RCTD](https://img.shields.io/badge/spacexr-v2.0-green)](https://github.com/dmcable/spacexr)
[![CARD](https://img.shields.io/badge/CARD-v1.1-orange)](https://github.com/SingleCellOpenData/CARD)
[![SPOTlight](https://img.shields.io/badge/SPOTlight-v1.6-purple)](https://bioconductor.org/packages/SPOTlight/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview
This repository provides an end-to-end bioinformatics pipeline for spatial transcriptomics analysis across multi-patient cohorts using **10x Visium** data. 

To overcome algorithm-specific biases and accurately resolve low-abundance / rare cell populations (such as Plasmablasts, Mast cells, or rare T-cell subsets), this workflow implements a **multi-engine deconvolution architecture**. It integrates multi-slice spatial dataset layers using **Seurat v5 (RPCA Integration)** and benchmarks single-cell RNA-sequencing (scRNA-seq) cell-type mapping across five statistical frameworks:
1. **RCTD (Full & Multi Mode):** Poisson / Overdispersion maximum-likelihood model.
2. **CARD:** Spatially informed Conditional Autoregressive (CAR) model incorporating $(x, y)$ coordinates.
3. **SPOTlight:** Seeded Non-Negative Matrix Factorization (NMF) with non-negative least squares (NNLS).
4. **SpatialDecon:** Constrained log-normal mixed-effects regression resistant to background noise.
5. **Seurat Anchor Transfer:** Canonical Correlation Analysis (CCA) probabilistic transfer.

## Key Features
* **Multi-Sample Cohort Integration:** Merges and normalizes 22+ Visium spatial slides (~50,000+ spatial spots; utilizes `SketchData()` for datasets scaling beyond 100k spots) using `SCTransform v2` and `JoinLayers` architecture in Seurat v5.
* **Multi-Algorithm Consensus Deconvolution:** Integrates R-native deconvolution methods to evaluate cell-type spatial distribution and produce robust cross-method consensus scores.
* **Rare Cell Type Sensitivity:** Employs **RCTD Full Mode** (unrestricted spot capacity) and **CARD Spatial Smoothing** to detect low-frequency populations without artificial truncation or false-negative masking.
* **In Situ Microenvironment Mapping:** Maps major and minor cell states (Cancer Epithelial, T-cells, Myeloid, CAFs, B-cells, Plasmablasts, PVL, Endothelial) directly onto tissue architecture without downsampling raw UMI counts.
* **Comparative Infiltration Analysis:** Evaluates spatial spatial niches and niche-specific cell density across demographic and clinical patient cohorts.

## Workflow Overview

```text
[Raw Visium Data] ──> [Quality Control] ──> [SCTransform Normalization]
                                                   │
                                                   ▼
[scRNA Reference] ──> [JoinLayers (Spatial)] ──> [Multi-Engine Deconvolution]
                                                   │
     ┌──────────────────┬──────────────────┬───────┴──────────┬──────────────────┐
     ▼                  ▼                  ▼                  ▼                  ▼
[RCTD Assay]      [CARD Assay]     [SPOTlight Assay] [SpatialDecon Assay] [SeuratTransfer Assay]
     │                  │                  │                  │                  │
     └──────────────────┴──────────────────┴───────┬──────────┴──────────────────┘
                                                   ▼
                                     [Consensus & Correlation]
                                                   │
                                                   ▼
                                      [Comparative Spatial Maps]
