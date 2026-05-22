# ============================================================
# 09_enrichment.R
# GO:BP + KEGG + GSEA enrichment analysis
# Yu G et al. (2012) doi:10.1089/omi.2011.0118
# Subramanian A et al. (2005) doi:10.1073/pnas.0506580102
# ============================================================

library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
library(ggplot2); library(writexl); library(dplyr)

hub_list  <- readRDS("data/hub_genes.rds")
all_hubs  <- unique(c(hub_list$lc, hub_list$as, hub_list$pn))
entrez    <- bitr(all_hubs, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# GO Biological Process
go_bp <- enrichGO(gene=entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP",
                   pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)
# KEGG
kegg  <- enrichKEGG(gene=entrez$ENTREZID, organism="hsa",
                     pAdjustMethod="BH", pvalueCutoff=0.05)

p_go   <- dotplot(go_bp, showCategory=15, title="GO BP — Hub Genes") + theme_bw(11)
p_kegg <- dotplot(kegg,  showCategory=15, title="KEGG — Hub Genes")  + theme_bw(11)
ggsave("results/figures/png/GO_BP_HubGenes.png",  p_go,   width=10, height=8, dpi=300)
ggsave("results/figures/png/KEGG_HubGenes.png",   p_kegg, width=10, height=8, dpi=300)
write_xlsx(as.data.frame(go_bp), "results/tables/GO_BP.xlsx")
write_xlsx(as.data.frame(kegg),  "results/tables/KEGG.xlsx")

# GSEA — Subramanian A et al. (2005)
deg_lc <- readRDS("data/deg_lc.rds")
run_gsea <- function(deg_df, disease_name) {
  ranked    <- deg_df %>% dplyr::filter(!is.na(logFC)) %>% arrange(desc(logFC))
  ranked    <- ranked[!duplicated(ranked$Gene), ]
  gene_list <- setNames(ranked$logFC, ranked$Gene)
  gsea_res  <- gseGO(geneList=gene_list, OrgDb=org.Hs.eg.db, keyType="SYMBOL",
                      ont="BP", minGSSize=10, maxGSSize=500,
                      pvalueCutoff=0.05, verbose=FALSE, eps=0)
  if (nrow(as.data.frame(gsea_res)) > 0) {
    p <- gseaplot2(gsea_res, geneSetID=1:min(3, nrow(as.data.frame(gsea_res))),
                   title=paste("GSEA —", disease_name))
    ggsave(paste0("results/figures/png/GSEA_", disease_name, ".png"),
           p, width=10, height=7, dpi=300)
  }
  return(gsea_res)
}
gsea_lc <- run_gsea(readRDS("data/deg_lc.rds"), "LungCancer")
gsea_as <- run_gsea(readRDS("data/deg_as.rds"), "Asthma")
gsea_pn <- run_gsea(readRDS("data/deg_pn.rds"), "Pneumonia")
cat("Step 09 complete. Enrichment analysis done.\n")
