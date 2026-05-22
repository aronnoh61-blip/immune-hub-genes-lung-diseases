# ============================================================
# 03_deg_analysis.R
# Differential expression analysis using limma
# Ritchie ME et al. (2015) doi:10.1093/nar/gkv007
# ============================================================

library(limma); library(ggplot2); library(ggrepel)

lc <- readRDS("data/lc_processed.rds")
as <- readRDS("data/as_processed.rds")
pn <- readRDS("data/pn_processed.rds")

# DEG function — adj.P < 0.05, |log2FC| > 0.5
run_deg <- function(expr_mat, group_vec, disease_name) {
  cat("\nRunning DEG:", disease_name, "\n")
  group  <- factor(group_vec, levels = c("Control", "Disease"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  fit  <- lmFit(expr_mat, design)
  cont <- makeContrasts(Disease - Control, levels = design)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2, trend = TRUE)
  deg  <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "B")
  deg$Gene   <- rownames(deg)
  deg$Status <- "Not Significant"
  deg$Status[deg$adj.P.Val < 0.05 & deg$logFC >  0.5] <- "Upregulated"
  deg$Status[deg$adj.P.Val < 0.05 & deg$logFC < -0.5] <- "Downregulated"
  cat("Up:", sum(deg$Status == "Upregulated"),
      "| Down:", sum(deg$Status == "Downregulated"), "\n")
  return(deg)
}

deg_lc <- run_deg(lc$expr, lc$group, "LungCancer")
deg_as <- run_deg(as$expr, as$group, "Asthma")
deg_pn <- run_deg(pn$expr, pn$group, "Pneumonia")

saveRDS(deg_lc, "data/deg_lc.rds")
saveRDS(deg_as, "data/deg_as.rds")
saveRDS(deg_pn, "data/deg_pn.rds")

# Volcano plots
make_volcano <- function(deg_df, disease_name) {
  top_genes <- head(deg_df[deg_df$Status != "Not Significant", ][
    order(deg_df[deg_df$Status != "Not Significant", "adj.P.Val"]), ], 15)
  colors <- c("Upregulated" = "#E24B4A", "Downregulated" = "#378ADD",
               "Not Significant" = "#B4B2A9")
  p <- ggplot(deg_df, aes(logFC, -log10(adj.P.Val), color = Status)) +
    geom_point(size = 1.5, alpha = 0.65) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray40") +
    geom_text_repel(data = top_genes, aes(label = Gene),
                    size = 3, color = "black", max.overlaps = 20) +
    scale_color_manual(values = colors) +
    labs(title = paste("Volcano Plot —", disease_name),
         x = "log2 Fold Change", y = "-log10 (adj. p-value)") +
    theme_bw(base_size = 13)
  ggsave(paste0("results/figures/png/Volcano_", disease_name, ".png"),
         p, width = 8, height = 7, dpi = 300)
  cat("Volcano plot saved:", disease_name, "\n")
  return(p)
}

make_volcano(deg_lc, "LungCancer")
make_volcano(deg_as, "Asthma")
make_volcano(deg_pn, "Pneumonia")

cat("Step 03 complete. DEG analysis done.\n")
