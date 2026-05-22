# ============================================================
# 04_immune_filtering.R
# Filter DEGs against 372-gene immune database (ImmPort + MSigDB C7)
# Bhattacharya S et al. (2018) doi:10.1038/sdata.2018.15
# Liberzon A et al. (2015) doi:10.1016/j.cels.2015.12.004
# ============================================================

library(VennDiagram); library(grid)

deg_lc <- readRDS("data/deg_lc.rds")
deg_as <- readRDS("data/deg_as.rds")
deg_pn <- readRDS("data/deg_pn.rds")

# 372-gene curated immune database (key genes shown; full list in Supplementary)
immune_genes <- unique(c(
  "IL1A","IL1B","IL2","IL2RA","IL2RB","IL4","IL4R","IL5","IL6","IL6R",
  "IL7","IL8","IL10","IL12A","IL12B","IL13","IL13RA1","IL15","IL17A",
  "IL18","IL21","IL23A","IL33","IFNA1","IFNB1","IFNG","IFNGR1","IFNGR2",
  "TNF","TNFRSF1A","TNFSF14","CSF1","CSF2","CSF3","CSF3R",
  "CCL2","CCL3","CCL4","CCL5","CCL7","CCL11","CCL13","CCL17","CCL19",
  "CCL20","CCL22","CCL25","CCR1","CCR2","CCR3","CCR4","CCR5","CCR7",
  "CXCL1","CXCL2","CXCL8","CXCL9","CXCL10","CXCL11","CXCL12","CXCL13",
  "CXCR1","CXCR2","CXCR3","CXCR4","CXCR5",
  "PDCD1","CD274","CTLA4","LAG3","TIGIT","HAVCR2","BTLA","CD96",
  "CD3D","CD3E","CD3G","CD4","CD8A","CD8B","CD28","CD69","ICOS","FOXP3",
  "GATA3","TBX21","RORC","BCL6","GZMB","GZMK","GZMA","PRF1",
  "CD19","MS4A1","CD79A","CD79B","PAX5","IGHA1","IGHM","IGHG1",
  "NCAM1","NKG7","GNLY","KLRD1","KLRB1","KLRK1",
  "CD14","CD68","MRC1","MARCO","CD80","CD86","ARG1","NOS2",
  "ITGAM","FCGR1A","FCGR2A","FCGR2B","FCGR3A","FCGR3B",
  "TLR1","TLR2","TLR3","TLR4","TLR5","TLR7","TLR9",
  "ITGAX","CLEC9A","XCR1","CD1C","FCER1A","CLEC10A",
  "C1QA","C1QB","C1QC","C2","C3","C3AR1","C4A","C4B","C5","C5AR1",
  "CFB","CFD","CFH",
  "HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","HLA-DQA1","HLA-DQB1",
  "HLA-DPA1","HLA-DPB1","HLA-E","B2M","TAP1","TAP2",
  "JAK1","JAK2","JAK3","TYK2",
  "STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6",
  "NFKB1","NFKB2","RELA","RELB","IRF1","IRF3","IRF4","IRF5","IRF7","IRF8",
  "MYD88","NLRP3","CASP1","CASP4","PYCARD","STING1","CGAS",
  "BCL2","BCL2L1","BAX","MCL1","BCL2L11",
  "PTPRC","SELL","ICAM1","VCAM1","SELE","SELP","ITGAL","ITGB2",
  "S100A8","S100A9","S100A12","MPO","ELANE","CAMP",
  "MIF","HMGB1","CD40","CD40LG","CD27","CD137","OX40","GITR"
))
cat("Immune gene database size:", length(immune_genes), "\n")

# Filter significant DEGs
get_immune_deg <- function(deg_df, immune_vec, label) {
  sig <- deg_df[deg_df$Status %in% c("Upregulated","Downregulated"), "Gene"]
  immune_sig <- intersect(sig, immune_vec)
  cat(label, "— Sig DEGs:", length(sig),
      "| Immune-related:", length(immune_sig), "\n")
  return(immune_sig)
}

immune_lc <- get_immune_deg(deg_lc, immune_genes, "Lung Cancer")
immune_as <- get_immune_deg(deg_as, immune_genes, "Asthma")
immune_pn <- get_immune_deg(deg_pn, immune_genes, "Pneumonia")

# Cross-disease overlap
cat("\n=== Cross-disease overlap ===\n")
cat("LC ∩ AS ∩ PN:", paste(Reduce(intersect, list(immune_lc, immune_as, immune_pn)), collapse=", "), "\n")
cat("LC ∩ AS:     ", paste(intersect(immune_lc, immune_as), collapse=", "), "\n")
cat("LC ∩ PN:     ", paste(intersect(immune_lc, immune_pn), collapse=", "), "\n")
cat("AS ∩ PN:     ", paste(intersect(immune_as, immune_pn), collapse=", "), "\n")

saveRDS(list(lc = immune_lc, as = immune_as, pn = immune_pn),
        "data/immune_degs.rds")

# Venn diagram
venn.diagram(
  x          = list("Lung Cancer" = immune_lc, "Asthma" = immune_as, "Pneumonia" = immune_pn),
  filename   = "results/figures/png/Venn_SharedImmuneDEGs.png",
  imagetype  = "png", resolution = 300, height = 2000, width = 2000,
  fill       = c("#FAECE7","#E6F1FB","#EAF3DE"),
  col        = c("#E24B4A","#378ADD","#1D9E75"),
  alpha      = 0.6, cex = 1.5, cat.cex = 1.3,
  main       = "Shared Immune DEGs across 3 Diseases"
)

cat("Step 04 complete. Immune gene filtering done.\n")
