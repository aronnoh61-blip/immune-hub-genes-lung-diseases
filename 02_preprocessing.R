# ============================================================
# 02_preprocessing.R
# Log2 transformation, quantile normalisation, gene mapping
# ============================================================

library(limma); library(GEOquery); library(sva)

# Gene mapping helper
map_genes <- function(expr_mat, feat_df, symbol_col) {
  feat_df <- feat_df[seq_len(nrow(expr_mat)), ]
  symbols <- gsub("\\s*///.*", "", feat_df[[symbol_col]])
  symbols <- trimws(symbols)
  keep    <- !is.na(symbols) & symbols != ""
  expr_mat <- expr_mat[keep, ]
  symbols  <- symbols[keep]
  rownames(expr_mat) <- symbols
  expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
  cat("Unique genes retained:", nrow(expr_mat), "\n")
  return(expr_mat)
}

# Log2 + quantile normalisation
preprocess <- function(expr_mat) {
  if (max(expr_mat, na.rm = TRUE) > 30) {
    expr_mat <- log2(expr_mat + 1)
    cat("Log2 transformation applied.\n")
  }
  expr_mat <- normalizeBetweenArrays(expr_mat, method = "quantile")
  cat("Quantile normalisation applied.\n")
  return(expr_mat)
}

# Lung Cancer (GSE19804)
gse_lc  <- readRDS("data/gse_lc_train.rds")
lc_expr <- preprocess(exprs(gse_lc))
lc_expr <- map_genes(lc_expr, fData(gse_lc), "Gene Symbol")
lc_group <- ifelse(grepl("tumor", pData(gse_lc)$source_name_ch1,
                          ignore.case = TRUE), "Disease", "Control")

# Asthma (GSE74986)
gse_as  <- readRDS("data/gse_as_train.rds")
as_expr <- preprocess(exprs(gse_as))
as_expr <- map_genes(as_expr, fData(gse_as), "GENE_SYMBOL")
as_group <- ifelse(grepl("HEALTHY", pData(gse_as)$source_name_ch1),
                   "Control", "Disease")

# Pneumonia (GSE65682) — bacterial + influenza A vs healthy
gse_pn   <- readRDS("data/gse_pn_train.rds")
pn_source <- pData(gse_pn)$source_name_ch1
keep_pn   <- grepl("bacterial pneumonia|influenza A pneumonia|healthy control",
                    pn_source, ignore.case = TRUE) &
             !grepl("mixed", pn_source, ignore.case = TRUE)
pn_expr  <- preprocess(exprs(gse_pn)[, keep_pn])
pn_expr  <- map_genes(pn_expr, fData(gse_pn), "Symbol")
pn_group <- ifelse(grepl("healthy control", pn_source[keep_pn],
                          ignore.case = TRUE), "Control", "Disease")

# Batch correction check (ComBat applied if needed after PCA)
# Johnson WE et al. (2007) doi:10.1093/biostatistics/kxj037

saveRDS(list(expr = lc_expr, group = lc_group), "data/lc_processed.rds")
saveRDS(list(expr = as_expr, group = as_group), "data/as_processed.rds")
saveRDS(list(expr = pn_expr, group = pn_group), "data/pn_processed.rds")

cat("Step 02 complete. Preprocessed data saved.\n")
cat("LC:", nrow(lc_expr), "genes,", ncol(lc_expr), "samples\n")
cat("AS:", nrow(as_expr), "genes,", ncol(as_expr), "samples\n")
cat("PN:", nrow(pn_expr), "genes,", ncol(pn_expr), "samples\n")
