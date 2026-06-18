# FL-MRS_Trajectory

**Integrative Computational Analysis of Public Fecal Lipidomics and Transcriptomics Datasets Suggests a Candidate Association Between the COX-2/CE(20:4) Axis and Colorectal Adenoma-Carcinoma Progression**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the complete computational pipeline for the manuscript:

> Tang B, Chen Y, Yang X, Gong H. *Integrative Computational Analysis of Public Fecal Lipidomics and Transcriptomics Datasets Suggests a Candidate Association Between the COX-2/CE(20:4) Axis and Colorectal Adenoma-Carcinoma Progression.* 

---

## Overview

We developed an **extreme-phenotype machine learning framework** to analyze publicly available fecal lipidomics data and identify candidate lipidomic features associated with colorectal adenoma-carcinoma progression. The pipeline integrates:

- **Fecal lipidomics** (Metabolomics Workbench: ST003798, ST002787)
- **Bulk transcriptomics** (TCGA-COAD)
- **Single-cell RNA-seq** (Broad Institute, c295)

Key outputs include a Fecal Lipidomic Malignancy Risk Score (FL-MRS), SHAP-based feature prioritization, cross-sectional pseudotime trajectory inference, and multiple sensitivity analyses.

**Important:** All findings are derived from publicly available retrospective data without independent experimental validation and should be regarded as hypothesis-generating.

---

## Repository Structure

| File | Description |
|------|-------------|
| `01_Main_Pipeline_ST003798.R` | Main analysis pipeline: data preprocessing, LASSO feature selection, Random Forest modeling, SHAP analysis, pseudotime inference, TCGA integration, single-cell analysis |
| `02_Supp_Basic_Metrics.R` | Supplementary: basic diagnostic performance metrics |
| `03_Supp_Advanced_Metrics_DCA.R` | Supplementary: calibration curves and decision curve analysis |
| `04_Supp_Feature_Extraction_and_AUC_CI.R` | Supplementary: feature extraction and bootstrap AUC confidence intervals |
| `05_Supp_Model_Comparisons_Table_S5.R` | Supplementary: multi-model comparison (RF vs SVM vs LR vs XGBoost) |
| `06_Supp_Final_Stats_and_TCGA_Covariates.R` | Supplementary: final summary statistics and TCGA covariate exploration |
| `07_Supp_Sensitivity_Analyses.R` | **New** – Supplementary: sensitivity analyses including OOB convergence, missing value imputation comparison (LOD/2), feature overlap between cohorts, pseudotime root node reversal, external validation pairwise comparisons, TCGA Wilcoxon validation, LASSO bootstrap stability analysis |

---

## Data Availability

All datasets used in this study are publicly available:

- **Fecal lipidomics (discovery)**: Metabolomics Workbench Study ID [ST003798](https://doi.org/10.21228/M8WR76)
- **Fecal lipidomics (external validation)**: Metabolomics Workbench Study ID [ST002787](https://doi.org/10.21228/M85X48)
- **Bulk transcriptomics**: TCGA-COAD via [GDC Data Portal](https://portal.gdc.cancer.gov/)
- **Single-cell RNA-seq**: Broad Institute Single Cell Portal, [c295](https://singlecell.broadinstitute.org/single_cell/study/SCP259)

---

## Reproducibility

All analyses were performed in **R version 4.3.1**. Required packages include:

`glmnet`, `randomForest`, `pROC`, `treeshap`, `princurve`, `Seurat` (v4.4.0), `ggplot2`, `VennDiagram` (optional)

To reproduce the full analysis:

1. Clone this repository
2. Download the required public datasets and place them in the appropriate directories (see individual scripts for paths)
3. Run `01_Main_Pipeline_ST003798.R` first to generate all core objects
4. Run supplementary scripts `02` through `07` in any order

**Note:** Script `07_Supp_Sensitivity_Analyses.R` generates supplementary figures S3, S4 and tables S6, S7.

---

## Supplementary Materials

The manuscript is accompanied by:

- **Supplementary Figures (PDF):** S1–S4 (calibration/DCA, study flowchart, OOB convergence, pseudotime root sensitivity)
- **Supplementary Tables (Excel):** S1–S11 (TCGA stats with Shapiro-Wilk, single-cell metrics, diagnostic performance, feature coefficients, model comparison, inter-cohort feature overlap, imputation sensitivity, Kruskal-Wallis tests, Dunn's post-hoc, Wilcoxon validation, bootstrap stability)

---

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

---

## Citation

If you use this code or the FL-MRS framework in your research, please cite:
Tang B, Chen Y, Yang X, Gong H. Integrative Computational Analysis of Public Fecal Lipidomics
and Transcriptomics Datasets Suggests a Candidate Association Between the COX-2/CE(20:4) Axis
and Colorectal Adenoma-Carcinoma Progression. 2026. GitHub: https://github.com/bingmoon/FL-MRS_Trajectory
