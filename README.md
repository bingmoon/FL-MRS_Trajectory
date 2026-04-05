# Machine Learning-Driven Fecal Lipidomic Trajectory (FL-MRS)

[![Journal](https://img.shields.io/badge/Submitted_to-Theranostics-red.svg)]()
[![Language](https://img.shields.io/badge/Language-R_4.3.1-blue.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)]()

This repository contains the official computational framework and reproducible analytical scripts for the manuscript:  
**"Machine learning-driven fecal lipidomic trajectory identifies the COX-2/CE(20:4) axis as a key early molecular correlate of colorectal adenoma-carcinoma transformation."**

## 📌 Project Overview
Accurately stratifying the malignant potential of histologically benign colorectal adenomas remains a major clinical challenge. To address this, we developed a novel **Extreme-Phenotype Machine Learning Framework** that establishes a pure diagnostic baseline by deliberately blinding the model to transitional adenoma stages during training. 

By projecting unseen adenomas onto this baseline and calculating the **Fecal Lipidomic Malignancy Risk Score (FL-MRS)**, we mathematically unmasked hidden heterogeneity, revealing that **53.4% of adenomas already harbor a high-risk, CRC-like lipidomic signature**. 

Integration of Explainable AI (TreeSHAP), pseudotime trajectory inference, and multi-omic validation (bulk & scRNA-seq) further pinpointed the **COX-2/CE(20:4) axis** as the primary driver of this early inflammatory burst within the Tumor Microenvironment (TME), providing a crucial temporal window for NSAID-based chemoprevention.

## 📂 Repository Structure & Reproducible Pipeline
The analysis is strictly modularized into 6 sequential R scripts to ensure 100% reproducibility. **Please execute them in numerical order.**

- `Scripts/`
  - **`01_Main_Pipeline_ST003798.R`**: The core analytical engine. Executes data preprocessing, LASSO-Random Forest extreme-phenotype training, SHAP interpretation, Pseudotime trajectory inference, and cross-cohort multi-omic (TCGA & scRNA-seq) target validation. (Generates Fig 1-9).
  - **`02_Supp_Basic_Metrics.R`**: Calculates comprehensive diagnostic metrics (Sensitivity, Specificity, PPV, NPV) for the independent test set.
  - **`03_Supp_Advanced_Metrics_DCA.R`**: Performs advanced clinical utility assessments, including Calibration Curves, Decision Curve Analysis (DCA), and 1,000x Bootstrap internal validation.
  - **`04_Supp_...R`**: *(Note: Please rename this description based on your actual Script 04, e.g., Supplementary visualizations or specific data processing).*
  - **`05_Supp_Model_Comparisons_Table_S5.R`**: Provides algorithmic justification by comparing the Extreme-Phenotype RF framework against "Mixed" models and other mainstream classifiers (XGBoost, SVM, Logistic Regression).
  - **`06_Supp_Final_Stats_and_TCGA_Covariates.R`**: Executes the ultimate statistical defense, including Hartigan's dip test for adenoma bimodality, external independent AUC calculation for CE(20:4), and Age/Sex-adjusted Multivariate Logistic Regression (GLM) for TCGA cohorts.

- `Data/`: Contains the core pre-processed fecal lipidomic matrices required to reproduce the machine learning framework.

## 💾 Data Availability
To facilitate immediate code reproducibility, the core lipidomic matrices are provided in the `Data/` directory:
- `MSdata_ST003798_1.txt` (Discovery Cohort)
- `st002787_positive_clean.txt` & `st002787_negative_clean.txt` (External Validation Cohort)

**Note on Massive Multi-Omic Data**: To comply with GitHub's file size limits (>100MB), the large-scale transcriptomic datasets are **not** hosted in this repository. They can be publicly accessed at:
1. **Tumor Microenvironment scRNA-seq**: [Broad Institute Single Cell Portal (c295)](https://singlecell.broadinstitute.org/single_cell/study/SCP1162/human-colon-cancer-atlas-c295).
2. **Host Transcriptomics**: TCGA-COAD bulk RNA-seq via the [GDC Data Portal](https://portal.gdc.cancer.gov/).

## 🛠️ Prerequisites
The scripts were built and optimized under **R version 4.3.1**. Ensure the following core packages are installed before running the pipeline:
```R
install.packages(c("glmnet", "randomForest", "xgboost", "e1071", "pROC", "caret", "dcurves", "diptest", "treeshap", "princurve"))
# For scRNA-seq spatial mapping:
install.packages("Seurat")
