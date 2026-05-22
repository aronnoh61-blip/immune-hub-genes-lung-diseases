# ============================================================
# 12_figures.R
# Generate all publication-ready figures
# ============================================================

library(pheatmap); library(ggplot2); library(RColorBrewer)
library(pROC); library(survival); library(survminer)

lc <- readRDS("data/lc_processed.rds")
as <- readRDS("data/as_processed.rds")
pn <- readRDS("data/pn_processed.rds")
hub_list <- readRDS("data/hub_genes.rds")

# Hub gene expression heatmaps
make_heatmap <- function(expr_mat, group_vec, hub_genes, disease_name) {
  hub_in  <- intersect(hub_genes, rownames(expr_mat))
  hub_expr <- expr_mat[hub_in, order(group_vec), drop=FALSE]
  ann_col  <- data.frame(Group=sort(group_vec))
  rownames(ann_col) <- colnames(hub_expr)
  png(paste0("results/figures/png/Heatmap_", disease_name, ".png"),
      width=3600, height=max(1200, length(hub_in)*300+600), res=300)
  pheatmap(hub_expr, annotation_col=ann_col, show_colnames=FALSE, scale="row",
           clustering_method="ward.D2",
           annotation_colors=list(Group=c(Disease="#E24B4A", Control="#378ADD")),
           color=colorRampPalette(c("#378ADD","white","#E24B4A"))(100),
           main=paste("Hub Gene Expression —", disease_name),
           fontsize_row=12, border_color=NA)
  dev.off()
  cat("Heatmap saved:", disease_name, "\n")
}
make_heatmap(lc$expr, lc$group, hub_list$lc, "LungCancer")
make_heatmap(as$expr, as$group, hub_list$as, "Asthma")

# AUC comparison bar chart
auc_df <- data.frame(
  Disease=c(rep("Lung Cancer",9), rep("Asthma",5), rep("Pneumonia",3)),
  Gene=c("IFNG","CCL13","C5","RELB","CXCL12","SELE","IGHA1","CD8A","HAVCR2",
         "HAVCR2","IL13RA1","MCL1","JAK2","CASP1","S100A12","IL2RB","CD96"),
  Training=c(0.930,0.941,0.922,0.935,0.927,0.891,0.876,0.860,0.942,
              0.917,0.893,0.937,0.910,0.940,0.994,0.999,0.998),
  Validation=c(0.601,0.740,0.702,0.585,0.873,0.784,0.615,0.767,0.747,
               0.507,0.567,0.708,0.637,0.560,1.000,0.998,1.000)
)
library(tidyr)
auc_long <- auc_df %>% pivot_longer(c(Training, Validation),
                                     names_to="Cohort", values_to="AUC")
p_auc <- ggplot(auc_long, aes(x=reorder(Gene,AUC), y=AUC, fill=Cohort)) +
  geom_bar(stat="identity", position=position_dodge(0.7), width=0.6) +
  geom_hline(yintercept=0.7, linetype="dashed", color="gray40") +
  coord_flip() + facet_wrap(~Disease, scales="free_y") +
  scale_fill_manual(values=c(Training="#378ADD", Validation="#E24B4A")) +
  scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1,0.2)) +
  labs(title="Hub Gene Diagnostic Performance — Training vs Validation",
       x=NULL, y="AUC", fill=NULL) +
  theme_bw(12) + theme(legend.position="bottom",
                        strip.text=element_text(size=12, face="bold"))
ggsave("results/figures/png/AUC_Comparison.png", p_auc, width=14, height=8, dpi=300)

cat("\nAll figures generated.\n")
cat("Step 12 complete.\n\n")
cat("================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("  Figures  →", normalizePath("results/figures/png"), "\n")
cat("  Tables   →", normalizePath("results/tables"), "\n")
cat("  Survival →", normalizePath("results/survival"), "\n")
cat("================================================\n")
