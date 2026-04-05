# Machine Learning-Driven Fecal Lipidomic Trajectory (FL-MRS)

[![Journal](https://img.shields.io/badge/Submitted_to-Theranostics-red.svg)]()
[![Language](https://img.shields.io/badge/Language-R_4.3.1-blue.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)]()

This repository contains the official computational framework and analytical scripts for the manuscript:  
**"Machine learning-driven fecal lipidomic trajectory identifies the COX-2/CE(20:4) axis as a key early molecular correlate of colorectal adenoma-carcinoma transformation."**

## 📌 Project Overview
Accurately stratifying the malignant potential of histologically benign colorectal adenomas remains a major clinical challenge. To address this, we developed a novel **Extreme-Phenotype Machine Learning Framework** that establishes a pure diagnostic baseline by deliberately blinding the model to transitional adenoma stages during training. 

By projecting unseen adenomas onto this baseline and calculating the **Fecal Lipidomic Malignancy Risk Score (FL-MRS)**, we mathematically unmasked hidden heterogeneity, revealing that **53.4% of adenomas already harbor a high-risk, CRC-like lipidomic signature**. 

Integration of Explainable AI (TreeSHAP), pseudotime trajectory inference, and multi-omic validation (bulk & scRNA-seq) further pinpointed the **COX-2/CE(20:4) axis** as the primary driver of this early inflammatory burst within the Tumor Microenvironment (TME), providing a crucial temporal window for NSAID-based chemoprevention.

## 📂 Repository Structure
- `Scripts/`: R scripts for the entire reproducible pipeline.
  - `01_Data_Preprocessing.R`: Data filtering, KNN imputation, and TIC normalization.
  - `02_Extreme_Phenotype_ML.R`: Implementation of the LASSO-Random Forest baseline excluding intermediate lesions.
  - `03_FL_MRS_Scoring_DCA.R`: Risk cutoff calculation, bootstrap internal validation, and Decision Curve Analysis (DCA).
  - `04_SHAP_Trajectory.R`: Non-linear trajectory inference and biomarker decryption via TreeSHAP.
  - `05_scRNAseq_Mapping.R`: Single-cell spatial mapping of the COX-2 (PTGS2) axis using Seurat.
- `Data/`: Contains the core pre-processed fecal lipidomic matrices required to reproduce the machine learning framework and trajectory analyses.

## 💾 Data Availability
To facilitate immediate code reproducibility, the core lipidomic matrices are provided in the `Data/` directory:
- `MSdata_ST003798_1.txt` (Discovery Cohort)
- `st002787_positive_clean.txt` & `st002787_negative_clean.txt` (External Validation Cohort)

**Note on Massive Multi-Omic Data**: To comply with GitHub's file size limits (>100MB), the large-scale transcriptomic datasets are **not** hosted in this repository. They can be publicly accessed at:
1. **Tumor Microenvironment scRNA-seq**: [Broad Institute Single Cell Portal (c295)](https://singlecell.broadinstitute.org/single_cell/study/SCP1162/human-colon-cancer-atlas-c295).
2. **Host Transcriptomics**: TCGA-COAD bulk RNA-seq via the [GDC Data Portal](https://portal.gdc.cancer.gov/).

## 🛠️ Prerequisites & Installation
The scripts were built and optimized under **R version 4.3.1**. Ensure the following core packages are installed before running the pipeline:
```R
install.packages(c("glmnet", "randomForest", "treeshap", "princurve", "ggplot2", "kknn"))
# For scRNA-seq spatial mapping:
install.packages("Seurat")
