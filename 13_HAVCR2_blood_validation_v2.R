# ============================================================
# Script 13 v2: HAVCR2 Blood-Based Asthma Validation
# Fixed: GPL loaded, correct label column used
# ============================================================

library(GEOquery); library(pROC); library(ggplot2)

cat("Loading GSE69683 with GPL annotation...\n")
gse <- getGEO("GSE69683", GSEMatrix = TRUE, getGPL = TRUE)
eset <- gse[[1]]

expr  <- exprs(eset)
pheno <- pData(eset)
fdata <- fData(eset)

# ── Check available label columns ────────────────────────────
cat("\nAvailable characteristic columns:\n")
char_cols <- grep("ch1", colnames(pheno), value=TRUE)
for (col in char_cols) {
  vals <- unique(pheno[[col]])
  if (length(vals) <= 10) cat(sprintf("  [%s]: %s\n", col, paste(vals, collapse=" | ")))
}

# ── Disease label: GSE69683 uses cohort:ch1 ──────────────────
cat("\nCohort values:\n")
print(table(pheno$`cohort:ch1`))

# asthma = GINA 1-4, healthy = GINA 0
pheno$group <- ifelse(grepl("GINA [1-4]|asthma", pheno$`cohort:ch1`,
                             ignore.case=TRUE), 1, 0)
cat("\nGroup distribution (1=asthma, 0=healthy):\n")
print(table(pheno$group))

# ── Find HAVCR2 probe ────────────────────────────────────────
cat("\nfData columns:", paste(colnames(fdata), collapse=", "), "\n")

# Search all columns
havcr2_probes <- c()
for (col in colnames(fdata)) {
  hits <- rownames(fdata)[grepl("^HAVCR2$", as.character(fdata[[col]]),
                                 ignore.case=TRUE)]
  if (length(hits) > 0) {
    cat(sprintf("Found in '%s': %s\n", col, paste(hits, collapse=", ")))
    havcr2_probes <- unique(c(havcr2_probes, hits))
  }
}

if (length(havcr2_probes) == 0) {
  # Last resort: search partial match in all columns
  cat("Trying partial match...\n")
  for (col in colnames(fdata)) {
    hits <- rownames(fdata)[grepl("HAVCR2", as.character(fdata[[col]]),
                                   ignore.case=TRUE)]
    if (length(hits) > 0) {
      cat(sprintf("Partial match in '%s': %s\n", col, paste(hits, collapse=", ")))
      havcr2_probes <- unique(c(havcr2_probes, hits))
    }
  }
}

if (length(havcr2_probes) == 0) stop("HAVCR2 not found even with partial match.")

# Best probe = highest mean expression
best_probe <- if (length(havcr2_probes) > 1) {
  means <- rowMeans(expr[havcr2_probes, , drop=FALSE])
  cat("Probe means:", paste(round(means,2), collapse=", "), "\n")
  names(which.max(means))
} else havcr2_probes

cat("Using probe:", best_probe, "\n")

# ── ROC Analysis ─────────────────────────────────────────────
havcr2_expr <- as.numeric(expr[best_probe, ])
roc_obj <- roc(pheno$group, havcr2_expr, levels=c(0,1),
               direction="<", ci=TRUE, boot.n=1000)
auc_val <- as.numeric(auc(roc_obj))
ci_vals <- as.numeric(ci(roc_obj))
wt      <- wilcox.test(havcr2_expr ~ pheno$group)

cat(sprintf("\n=== RESULTS ===\nAUC = %.3f (95%% CI: %.3f-%.3f)\np = %.4f\n",
            auc_val, ci_vals[1], ci_vals[3], wt$p.value))

# ── ROC Plot ─────────────────────────────────────────────────
roc_df <- data.frame(FPR=1-roc_obj$specificities, TPR=roc_obj$sensitivities)
p <- ggplot(roc_df, aes(x=FPR, y=TPR)) +
  geom_line(color="#8E44AD", linewidth=1.8) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="#BDC3C7") +
  annotate("text", x=0.6, y=0.15,
           label=sprintf("AUC = %.3f\n95%% CI: %.3f-%.3f\np = %.4f",
                         auc_val, ci_vals[1], ci_vals[3], wt$p.value),
           size=4.5, hjust=0) +
  labs(title="HAVCR2 — Blood-Based Asthma Validation",
       subtitle="GSE69683: Whole Blood | Asthma vs Healthy",
       x="1 - Specificity", y="Sensitivity") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold")) + coord_equal()

ggsave("HAVCR2_blood_validation_ROC.png", p, width=6, height=6, dpi=300, bg="white")

# ── Boxplot ──────────────────────────────────────────────────
box_df <- data.frame(Expression=havcr2_expr,
                     Group=ifelse(pheno$group==1,"Asthma","Healthy"))
p2 <- ggplot(box_df, aes(x=Group, y=Expression, fill=Group)) +
  geom_boxplot(width=0.5, alpha=0.8, outlier.shape=21) +
  geom_jitter(width=0.12, size=0.8, alpha=0.3) +
  scale_fill_manual(values=c("Asthma"="#8E44AD","Healthy"="#27AE60")) +
  labs(title="HAVCR2 Expression — GSE69683",
       subtitle=sprintf("Wilcoxon p = %.4f", wt$p.value),
       x=NULL, y="Expression (log2)") +
  theme_classic(base_size=13) + theme(legend.position="none",
                                       plot.title=element_text(face="bold"))
ggsave("HAVCR2_blood_validation_boxplot.png", p2, width=5, height=5, dpi=300, bg="white")

# ── Save CSV ─────────────────────────────────────────────────
write.csv(data.frame(Dataset="GSE69683", Tissue="Whole Blood", Gene="HAVCR2",
                     N_asthma=sum(pheno$group==1), N_control=sum(pheno$group==0),
                     AUC=round(auc_val,3), CI_lower=round(ci_vals[1],3),
                     CI_upper=round(ci_vals[3],3), Wilcoxon_p=round(wt$p.value,4)),
          "HAVCR2_blood_validation_results.csv", row.names=FALSE)

cat("\n✅ Done! Send these 3 files:\n")
cat("  1. HAVCR2_blood_validation_ROC.png\n")
cat("  2. HAVCR2_blood_validation_boxplot.png\n")
cat("  3. HAVCR2_blood_validation_results.csv\n")
