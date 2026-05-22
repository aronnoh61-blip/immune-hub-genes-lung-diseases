# ============================================================
# 01_data_download.R
# Download GEO datasets (Training + Validation cohorts)
# ============================================================

library(GEOquery)
options(timeout = 600)

cat("Downloading training datasets...\n")

# Training cohorts
gse_lc <- getGEO("GSE19804", GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = "data/")[[1]]
gse_as <- getGEO("GSE74986", GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = "data/")[[1]]
gse_pn <- getGEO("GSE65682", GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = "data/")[[1]]

saveRDS(gse_lc, "data/gse_lc_train.rds")
saveRDS(gse_as, "data/gse_as_train.rds")
saveRDS(gse_pn, "data/gse_pn_train.rds")
cat("Training datasets saved.\n")

# Validation cohorts
cat("Downloading validation datasets...\n")
gse_lc_val <- getGEO("GSE30219", GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = "data/")[[1]]
gse_as_val <- getGEO("GSE41861", GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = "data/")[[1]]
gse_pn_val <- getGEO("GSE40012", GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = "data/")[[1]]

saveRDS(gse_lc_val, "data/gse_lc_val.rds")
saveRDS(gse_as_val, "data/gse_as_val.rds")
saveRDS(gse_pn_val, "data/gse_pn_val.rds")
cat("Validation datasets saved.\n")

# TCGA-LUSC expression data
cat("Downloading TCGA-LUSC expression data...\n")
url_lusc <- paste0("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/",
                   "download/TCGA.LUSC.sampleMap%2FHiSeqV2.gz")
download.file(url_lusc, destfile = "data/TCGA.LUSC.HiSeqV2.gz", mode = "wb")

url_luad <- paste0("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/",
                   "download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz")
download.file(url_luad, destfile = "data/TCGA.LUAD.HiSeqV2.gz", mode = "wb")
cat("TCGA data downloaded.\n")

cat("Step 01 complete. All datasets downloaded.\n")
