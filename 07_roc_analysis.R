# ============================================================
# 07_roc_analysis.R
# Diagnostic ROC analysis with 1000-bootstrap CIs
# Robin X et al. (2011) doi:10.1186/1471-2105-12-77
# ============================================================

library(pROC); library(ggplot2); library(writexl)

lc <- readRDS("data/lc_processed.rds")
as <- readRDS("data/as_processed.rds")
pn <- readRDS("data/pn_processed.rds")
hub_list <- readRDS("data/hub_genes.rds")

run_roc <- function(expr_mat, group_vec, hub_genes, disease_name) {
  cat("\nROC:", disease_name, "\n")
  genes_in <- intersect(hub_genes, rownames(expr_mat))
  y_real   <- as.numeric(group_vec == "Disease")
  roc_list <- list()
  for (gene in genes_in) {
    roc_obj <- roc(y_real, as.numeric(expr_mat[gene, ]), quiet=TRUE,
                   ci=TRUE, ci.method="bootstrap", boot.n=1000)
    roc_list[[gene]] <- data.frame(
      Gene=gene, AUC=round(as.numeric(auc(roc_obj)), 3),
      CI_low=round(as.numeric(ci(roc_obj))[1], 3),
      CI_high=round(as.numeric(ci(roc_obj))[3], 3))
    cat(gene, "AUC:", round(as.numeric(auc(roc_obj)), 3), "\n")
  }
  roc_df <- do.call(rbind, roc_list)
  write_xlsx(roc_df, paste0("results/tables/ROC_", disease_name, ".xlsx"))
  return(roc_df)
}

roc_lc <- run_roc(lc$expr, lc$group, hub_list$lc, "LungCancer")
roc_as <- run_roc(as$expr, as$group, hub_list$as, "Asthma")
roc_pn <- run_roc(pn$expr, pn$group, hub_list$pn, "Pneumonia")

# External validation ROC
gse_lc_val  <- readRDS("data/gse_lc_val.rds")
gse_as_val  <- readRDS("data/gse_as_val.rds")
gse_pn_val  <- readRDS("data/gse_pn_val.rds")

source("R/02_preprocessing.R")   # load map_genes and preprocess helpers

lc_val_expr  <- preprocess(exprs(gse_lc_val))
lc_val_expr  <- map_genes(lc_val_expr, fData(gse_lc_val), "Gene Symbol")
lc_val_group <- ifelse(grepl("Non Tumoral", pData(gse_lc_val)$source_name_ch1,
                              ignore.case=TRUE), "Control", "Disease")

as_val_expr  <- preprocess(exprs(gse_as_val))
as_val_expr  <- map_genes(as_val_expr, fData(gse_as_val), "Gene Symbol")
as_val_group <- ifelse(grepl("Healthy Control", pData(gse_as_val)$source_name_ch1,
                              ignore.case=TRUE), "Control", "Disease")

pn_val_expr  <- preprocess(exprs(gse_pn_val))
pn_val_expr  <- map_genes(pn_val_expr, fData(gse_pn_val), "Symbol")
pn_keep      <- grepl("bacterial pneumonia.*Day_1|influenza.*Day_1|healthy control.*Day_1",
                       pData(gse_pn_val)$source_name_ch1, ignore.case=TRUE)
pn_val_expr  <- pn_val_expr[, pn_keep]
pn_val_group <- ifelse(grepl("healthy control", pData(gse_pn_val)$source_name_ch1[pn_keep],
                              ignore.case=TRUE), "Control", "Disease")

roc_lc_val <- run_roc(lc_val_expr, lc_val_group, hub_list$lc, "LungCancer_Validation")
roc_as_val <- run_roc(as_val_expr, as_val_group, hub_list$as, "Asthma_Validation")
roc_pn_val <- run_roc(pn_val_expr, pn_val_group, hub_list$pn, "Pneumonia_Validation")

cat("Step 07 complete. ROC analysis done.\n")
