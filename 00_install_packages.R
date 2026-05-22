# ============================================================
# 00_install_packages.R
# Install all required packages
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# CRAN packages
install.packages(c(
  "tidyverse", "ggplot2", "pheatmap", "RColorBrewer",
  "VennDiagram", "glmnet", "randomForest", "e1071",
  "pROC", "caret", "ggpubr", "cowplot", "igraph",
  "ggraph", "writexl", "openxlsx", "corrplot",
  "ggrepel", "scales", "reshape2", "gridExtra", "pdftools"
), dependencies = TRUE)

# Bioconductor packages
BiocManager::install(c(
  "GEOquery", "limma", "edgeR", "DESeq2", "sva",
  "clusterProfiler", "org.Hs.eg.db", "DOSE",
  "enrichplot", "STRINGdb", "GSEABase", "GSVA",
  "survival", "survminer", "TCGAbiolinks",
  "AnnotationDbi", "biomaRt"
), ask = FALSE, update = FALSE)

cat("All packages installed successfully.\n")

# Load all libraries to verify
suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(sva)
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
  library(STRINGdb); library(GSVA); library(survival); library(survminer)
  library(tidyverse); library(ggplot2); library(pheatmap)
  library(RColorBrewer); library(VennDiagram); library(glmnet)
  library(randomForest); library(e1071); library(pROC); library(caret)
  library(igraph); library(ggraph); library(writexl); library(corrplot)
  library(ggrepel); library(reshape2)
})

cat("All libraries loaded successfully.\n")

# Create output folders
folders <- c("data", "results/figures", "results/figures/png",
             "results/tables", "results/survival")
sapply(folders, dir.create, showWarnings = FALSE, recursive = TRUE)
cat("Output folders created.\n")
