# Shared Immune Biomarkers Across Malignant, Allergic, and Infectious Lung Diseases

[![R Version](https://img.shields.io/badge/R-4.3.2-blue)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)

## Overview

This repository contains all R analysis scripts for the manuscript:

> **"Shared Immune Biomarkers Across Malignant, Allergic, and Infectious Lung Diseases: An Integrative Bioinformatics and Machine Learning Framework"**
>
> Md. Jamil Hassan, Md. Kamal Uddin
> Department of Statistics and Data Science, Islamic University, Kushtia-7003, Bangladesh

---

## Key Findings

| Disease | Hub Genes | Best AUC |
|---------|-----------|----------|
| Lung Cancer | IFNG, CCL13, C5, RELB, CXCL12, SELE, IGHA1, CD8A, **HAVCR2** | 0.942 (HAVCR2) |
| Asthma | **HAVCR2**, IL13RA1, MCL1, JAK2, CASP1 | 0.940 (CASP1) |
| Pneumonia | S100A12, IL2RB, CD96 | 0.999 (IL2RB) |

- **HAVCR2 (TIM-3)** identified as shared cross-disease immune checkpoint hub (Lung Cancer ∩ Asthma)
- **CXCL12** associated with improved overall survival in TCGA-LUSC (HR = 0.74, p = 0.032)
- **JAK–STAT signalling** is the most significantly enriched pathway (adj. p = 0.0016)

---

## Repository Structure

```
lung-immune-biomarkers/
├── README.md
├── LICENSE
├── R/
│   ├── 00_install_packages.R     # Package installation
│   ├── 01_data_download.R        # GEO dataset download
│   ├── 02_preprocessing.R        # Normalization, gene mapping
│   ├── 03_deg_analysis.R         # Differential expression (limma)
│   ├── 04_immune_filtering.R     # Immune gene filtering + Venn
│   ├── 05_ppi_network.R          # STRING PPI network construction
│   ├── 06_machine_learning.R     # LASSO + Random Forest + SVM-RFE
│   ├── 07_roc_analysis.R         # Diagnostic performance (pROC)
│   ├── 08_survival_analysis.R    # TCGA survival (LUAD + LUSC)
│   ├── 09_enrichment.R           # GO/KEGG + GSEA enrichment
│   ├── 10_immune_infiltration.R  # ssGSEA immune cell scoring
│   ├── 11_drug_target.R          # Drug–target interaction analysis
│   └── 12_figures.R              # Publication-ready figure generation
└── data/                         # Created automatically at runtime
    └── (GEO datasets downloaded here)
```

---

## Datasets

### Training Cohorts (NCBI GEO)
| Disease | GEO ID | Platform | Samples |
|---------|--------|----------|---------|
| Lung Cancer | GSE19804 | Affymetrix HG-U133 Plus 2.0 | 60 tumour + 60 normal |
| Asthma | GSE74986 | Agilent | 74 asthmatic + 12 healthy |
| Pneumonia | GSE65682 | Illumina HumanHT-12 v4 | 100 pneumonia + 36 healthy |

### External Validation Cohorts
| Disease | GEO ID | Samples |
|---------|--------|---------|
| Lung Cancer | GSE30219 | n = 307 |
| Asthma | GSE41861 | n = 138 |
| Pneumonia | GSE40012 | n = 42 |

### Survival Data
- **TCGA-LUAD** (n = 503) and **TCGA-LUSC** (n = 494) from [UCSC Xena](https://xenabrowser.net)

---

## Requirements

### R Version
```
R >= 4.3.2
```

### Key Packages

**CRAN:**
```r
tidyverse, ggplot2, pheatmap, RColorBrewer, VennDiagram,
glmnet, randomForest, e1071, pROC, caret, ggpubr, cowplot,
igraph, ggraph, writexl, corrplot, ggrepel, reshape2
```

**Bioconductor:**
```r
GEOquery, limma, sva, clusterProfiler, org.Hs.eg.db,
enrichplot, STRINGdb, GSVA, survival, survminer, TCGAbiolinks
```

Install all packages by running:
```r
source("R/00_install_packages.R")
```

---

## How to Run

Run scripts in order (01 → 12):

```r
# Step 1: Install packages
source("R/00_install_packages.R")

# Step 2: Download GEO data
source("R/01_data_download.R")

# Step 3: Preprocess
source("R/02_preprocessing.R")

# Step 4-12: Run remaining analyses
# See individual scripts for details
```

**Expected runtime:** ~3–5 hours (data download dependent on internet speed)

---

## Output Files

```
results/
├── figures/          # Volcano plots, PPI networks, ROC curves,
│   │                 # Heatmaps, KM survival plots
│   └── png/          # PNG versions (300 DPI) for publication
├── tables/           # Hub gene summary, ROC tables, survival results,
│                     # Drug-target interactions, enrichment results
└── survival/         # Kaplan-Meier plots per gene
```

---

## Citation

If you use these scripts, please cite:

> Hassan MJ, Uddin MK. Shared Immune Biomarkers Across Malignant, Allergic, and Infectious Lung Diseases: An Integrative Bioinformatics and Machine Learning Framework. *[Journal Name]*, 2025.

---

## Contact

**Corresponding author:** Md. Kamal Uddin
**Email:** kamal@iu.ac.bd
**Institution:** Department of Statistics and Data Science, Islamic University, Kushtia-7003, Bangladesh

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
