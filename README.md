# Spatial-Project
# Spatial Transcriptomics Architecture & Deconvolution Pipeline

[![Seurat v5](https://img.shields.io/badge/Seurat-v5.0-blue)](https://satijalab.org/seurat/)
[![RCTD / spacexr](https://img.shields.io/badge/spacexr-v2.0-green)](https://github.com/dmcable/spacexr)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview
This repository provides an end-to-end bioinformatics pipeline for spatial transcriptomics analysis across multi-patient cohorts using **10x Visium** data. 

The workflow integrates multi-slice spatial dataset layers using **Seurat v5 (RPCA Integration)** and performs cell-type deconvolution via **Robust Cell Type Deconvolution (RCTD / spacexr)** using an un-normalized single-cell RNA-sequencing (scRNA-seq) reference.

## Key Features
* **Multi-Sample Integration:** Merges and normalizes 22 + Visium spatial slides (~50000+ spatial spots: If more than 100k spots use sketch assay using SketchData() function) using `SCTransform v2` and `JoinLayers` architecture in Seurat v5.
* **Cell-Type Mapping In Situ:** Maps high-resolution cell types (Cancer Epithelial, T-cells, Myeloid, CAFs, B-cells) to spatial tissue coordinates without downsampling raw UMI counts.
* **Comparative Microenvironment Profiling:** Evaluates spatial infiltration patterns across demographic and clinical patient cohorts.

## Workflow Overview

```text
[Raw Visium Data] ──> [Quality Control] ──> [SCTransform Normalization]
                                                    │
                                                    ▼
[scRNA Reference] ──> [RCTD Deconvolution] <─── [JoinLayers (Spatial)]
                               │
                               ▼
                    [RCTD Proportions Assay] ──> [Comparative Spatial Plots]
