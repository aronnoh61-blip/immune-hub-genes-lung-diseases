# ============================================================
# 06_machine_learning.R
# LASSO + Random Forest + SVM-RFE consensus hub gene selection
# Tibshirani R (1996) doi:10.1111/j.2517-6161.1996.tb02080.x
# Breiman L (2001) doi:10.1023/A:1010933404324
# Guyon I et al. (2002) doi:10.1023/A:1012487302797
# ============================================================

library(glmnet); library(randomForest); library(caret); library(e1071); library(dplyr)

lc <- readRDS("data/lc_processed.rds")
as <- readRDS("data/as_processed.rds")
pn <- readRDS("data/pn_processed.rds")
immune_degs <- readRDS("data/immune_degs.rds")

# Use union of all immune DEGs as input pool
all_immune <- unique(c(immune_degs$lc, immune_degs$as, immune_degs$pn))

run_lasso <- function(expr_mat, group_vec, gene_vec, disease_name) {
  cat("\nLASSO:", disease_name, "\n")
  genes_in <- intersect(gene_vec, rownames(expr_mat))
  X <- t(expr_mat[genes_in, ])
  y <- as.numeric(group_vec == "Disease")
  set.seed(123)
  cv_fit  <- cv.glmnet(X, y, alpha=1, family="binomial", nfolds=10, type.measure="auc")
  coefs   <- coef(cv_fit, s="lambda.min")
  coef_df <- data.frame(Gene=rownames(coefs), Coef=as.numeric(coefs)) %>%
    dplyr::filter(Gene != "(Intercept)" & Coef != 0)
  selected <- intersect(coef_df$Gene, gene_vec)
  cat("Selected:", length(selected), "—", paste(selected, collapse=", "), "\n")
  return(selected)
}

run_rf <- function(expr_mat, group_vec, gene_vec, disease_name) {
  cat("\nRandom Forest:", disease_name, "\n")
  genes_in <- intersect(gene_vec, rownames(expr_mat))
  X <- data.frame(t(expr_mat[genes_in, ]))
  y <- factor(ifelse(group_vec == "Disease", "Disease", "Control"))
  set.seed(123)
  rf_fit <- randomForest(x=X, y=y, ntree=1000, mtry=floor(sqrt(ncol(X))), importance=TRUE)
  imp_df  <- data.frame(Gene=rownames(importance(rf_fit)),
                         MDA=importance(rf_fit)[,"MeanDecreaseAccuracy"]) %>%
    arrange(desc(MDA))
  top_genes <- head(imp_df$Gene, 10)
  cat("Top 10:", paste(top_genes, collapse=", "), "\n")
  return(top_genes)
}

run_svm_rfe <- function(expr_mat, group_vec, gene_vec, disease_name) {
  cat("\nSVM-RFE:", disease_name, "\n")
  genes_in <- intersect(gene_vec, rownames(expr_mat))
  X <- data.frame(t(expr_mat[genes_in, ]))
  y <- factor(ifelse(group_vec == "Disease", "Disease", "Control"))
  ctrl <- rfeControl(functions=caretFuncs, method="cv", number=5, verbose=FALSE)
  sizes <- c(3,5,8,10)[c(3,5,8,10) <= ncol(X)]
  set.seed(123)
  rfe_res   <- rfe(x=X, y=y, sizes=sizes, rfeControl=ctrl, method="svmLinear")
  svm_genes <- predictors(rfe_res)
  cat("Selected:", length(svm_genes), "—", paste(svm_genes, collapse=", "), "\n")
  return(svm_genes)
}

lasso_lc <- run_lasso(lc$expr, lc$group, all_immune, "LungCancer")
lasso_as <- run_lasso(as$expr, as$group, all_immune, "Asthma")
lasso_pn <- run_lasso(pn$expr, pn$group, all_immune, "Pneumonia")

rf_lc <- run_rf(lc$expr, lc$group, all_immune, "LungCancer")
rf_as <- run_rf(as$expr, as$group, all_immune, "Asthma")
rf_pn <- run_rf(pn$expr, pn$group, all_immune, "Pneumonia")

svm_lc <- run_svm_rfe(lc$expr, lc$group, all_immune, "LungCancer")
svm_as <- run_svm_rfe(as$expr, as$group, all_immune, "Asthma")
svm_pn <- run_svm_rfe(pn$expr, pn$group, all_immune, "Pneumonia")

# Consensus: selected by >= 2 of 3 algorithms
get_hub <- function(l, r, s) unique(c(intersect(l,r), intersect(l,s), intersect(r,s)))
hub_lc <- get_hub(lasso_lc, rf_lc, svm_lc)
hub_as <- get_hub(lasso_as, rf_as, svm_as)
hub_pn <- get_hub(lasso_pn, rf_pn, svm_pn)

cat("\n=== CONSENSUS HUB GENES ===\n")
cat("Lung Cancer:", paste(hub_lc, collapse=", "), "\n")
cat("Asthma:     ", paste(hub_as, collapse=", "), "\n")
cat("Pneumonia:  ", paste(hub_pn, collapse=", "), "\n")
cat("Shared LC∩AS:", paste(intersect(hub_lc, hub_as), collapse=", "), "\n")

saveRDS(list(lc=hub_lc, as=hub_as, pn=hub_pn), "data/hub_genes.rds")
cat("Step 06 complete. Hub genes identified.\n")
