# ============================================================
# 10_immune_infiltration.R
# ssGSEA immune cell infiltration via GSVA
# Barbie DA et al. (2009) doi:10.1038/nature08460
# Hanzelmann S et al. (2013) doi:10.1186/1471-2105-14-7
# ============================================================

library(GSVA); library(pheatmap); library(RColorBrewer); library(writexl)

lc <- readRDS("data/lc_processed.rds")
as <- readRDS("data/as_processed.rds")
pn <- readRDS("data/pn_processed.rds")
hub_list <- readRDS("data/hub_genes.rds")

# 10 immune cell gene sets (simplified; full sets from MSigDB C7)
immune_sets <- list(
  CD8_T_cells  = c("CD8A","CD8B","CD3D","CD3E","GZMB","PRF1"),
  CD4_T_cells  = c("CD4","CD3D","CD3E","IL2","ICOS"),
  Treg         = c("FOXP3","IL2RA","CTLA4","IKZF2"),
  NK_cells     = c("NCAM1","NKG7","GNLY","KLRD1","KLRK1"),
  B_cells      = c("CD19","MS4A1","CD79A","CD79B","PAX5"),
  Macrophages  = c("CD68","MRC1","MARCO","CD80","CD86"),
  Dendritic    = c("ITGAX","CLEC9A","CD1C","FCER1A"),
  Neutrophils  = c("FCGR3B","S100A8","S100A9","MPO","ELANE"),
  Mast_cells   = c("KIT","MS4A2","FCER1A","CPA3"),
  Monocytes    = c("CD14","FCGR2A","TLR2","TLR4","CSF1R")
)

run_ssgsea <- function(expr_mat, group_vec, disease_name) {
  param <- gsvaParam(expr_mat, immune_sets, method="ssgsea")
  scores <- gsva(param, verbose=FALSE)
  diff_scores <- scores[, order(group_vec)]
  ann_col    <- data.frame(Group=sort(group_vec))
  rownames(ann_col) <- colnames(diff_scores)
  pdf(paste0("results/figures/png/ssGSEA_", disease_name, ".png"),
      width=14, height=5)
  pheatmap(diff_scores, annotation_col=ann_col, show_colnames=FALSE,
           color=colorRampPalette(c("#378ADD","white","#E24B4A"))(100),
           main=paste("ssGSEA Immune Infiltration —", disease_name),
           border_color=NA, fontsize_row=11)
  dev.off()
  write_xlsx(as.data.frame(scores), paste0("results/tables/ssGSEA_", disease_name, ".xlsx"))
  return(scores)
}

imm_lc <- run_ssgsea(lc$expr, lc$group, "LungCancer")
imm_as <- run_ssgsea(as$expr, as$group, "Asthma")
imm_pn <- run_ssgsea(pn$expr, pn$group, "Pneumonia")
cat("Step 10 complete. Immune infiltration done.\n")
