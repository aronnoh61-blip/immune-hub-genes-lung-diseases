# ============================================================
# 05_ppi_network.R
# PPI network via STRINGdb v11.5 (score >= 400)
# Szklarczyk D et al. (2021) doi:10.1093/nar/gkaa1074
# igraph: Csardi G & Nepusz T (2006) InterJournal 1695
# ============================================================

library(STRINGdb); library(igraph); library(ggraph); library(ggplot2); library(dplyr)

string_db <- STRINGdb$new(version="11.5", species=9606,
                           score_threshold=400, input_directory="data/")
immune_degs <- readRDS("data/immune_degs.rds")

build_ppi <- function(gene_vec, disease_name) {
  cat("\nBuilding PPI:", disease_name, "\n")
  gene_df <- data.frame(Gene = gene_vec)
  mapped  <- string_db$map(gene_df, "Gene", removeUnmappedRows = TRUE)
  cat("Mapped:", nrow(mapped), "/", length(gene_vec), "\n")
  if (nrow(mapped) < 3) return(NULL)
  interactions <- string_db$get_interactions(mapped$STRING_id)
  if (nrow(interactions) == 0) return(NULL)
  id2gene <- setNames(mapped$Gene, mapped$STRING_id)
  interactions$from_gene <- id2gene[interactions$from]
  interactions$to_gene   <- id2gene[interactions$to]
  interactions <- interactions[!is.na(interactions$from_gene) &
                                 !is.na(interactions$to_gene), ]
  g <- simplify(graph_from_data_frame(
    interactions[, c("from_gene","to_gene","combined_score")], directed=FALSE))
  V(g)$degree      <- degree(g)
  V(g)$betweenness <- betweenness(g, normalized=TRUE)
  metrics <- data.frame(Gene=V(g)$name, Degree=V(g)$degree,
                         Betweenness=round(V(g)$betweenness, 4)) %>%
    arrange(desc(Degree))
  cat("Nodes:", vcount(g), "| Edges:", ecount(g), "\n")
  cat("Top 5:", paste(head(metrics$Gene, 5), collapse=", "), "\n")
  set.seed(42)
  p <- ggraph(g, layout="fr") +
    geom_edge_link(color="gray70", alpha=0.6, width=0.4) +
    geom_node_point(aes(size=degree, color=degree)) +
    geom_node_text(aes(label=name), size=3, repel=TRUE) +
    scale_color_gradient(low="#B5D4F4", high="#E24B4A") +
    scale_size_continuous(range=c(2,10)) + theme_void() +
    labs(title=paste("PPI Network —", disease_name), color="Degree", size="Degree")
  ggsave(paste0("results/figures/png/PPI_", disease_name, ".png"),
         p, width=10, height=8, dpi=300)
  return(list(graph=g, metrics=metrics))
}

ppi_lc <- build_ppi(immune_degs$lc, "LungCancer")
ppi_as <- build_ppi(immune_degs$as, "Asthma")
ppi_pn <- build_ppi(immune_degs$pn, "Pneumonia")

saveRDS(list(lc=ppi_lc, as=ppi_as, pn=ppi_pn), "data/ppi_results.rds")

candidates <- list(lc=head(ppi_lc$metrics$Gene, 20),
                   as=head(ppi_as$metrics$Gene, 20),
                   pn=head(ppi_pn$metrics$Gene, 20))
saveRDS(candidates, "data/candidates.rds")
cat("Step 05 complete. PPI networks built.\n")
