#!/usr/bin/env Rscript
# 07_Fig2D_enrichr_network.R
# Figure 2D: miR-4532 target gene enrichment via enrichR
# If FDR<0.1 terms exist: network plot via plot_mirna_network()
# If no terms survive FDR<0.1: fallback dotplot of top nominal terms
# Output: results/Fig2D_enrichr_network.pdf OR results/Fig2D_enrichr_dotplot.png
#
# Usage:
#   Rscript src/07_Fig2D_enrichr_network.R

library(here)
library(tidyverse)
library(DBI)
library(RSQLite)
library(enrichR)

source(here("src", "utils.R"))

# --- Query mir2gene.sqlite for hsa-mir-4532 targets ---
cat("Querying mir2gene.sqlite for hsa-mir-4532 targets...\n")
mydb <- dbConnect(RSQLite::SQLite(), here("data", "annotations", "mir2gene.sqlite"))
mir2gene <- dbGetQuery(mydb, 'SELECT * FROM hsa')
dbDisconnect(mydb)

mirdat <- mir2gene %>% dplyr::filter(mir_id %in% "hsa-mir-4532")
cat("Number of target genes for hsa-mir-4532:", nrow(mirdat), "\n")

# --- Run enrichR ---
dbs <- c("GO_Biological_Process_2025",
         "GO_Cellular_Component_2025",
         "GO_Molecular_Function_2025",
         "WikiPathways_2024_Human")

cat("Running enrichR with", length(dbs), "databases...\n")
result <- enrichR::enrichr(mirdat$symbol, dbs)

# --- Filter for FDR < 0.1 ---
go <- lapply(c("hsa-mir-4532"), function(mir) {
  mirdat <- mir2gene %>% filter(mir_id %in% mir)
  result <- enrichR::enrichr(mirdat$symbol, dbs)
  df <- lapply(names(result), function(i) {
    subset(result[[i]], Adjusted.P.value < 0.1) %>%
      mutate(db = i)
  }) %>%
    do.call(rbind, .)
  df$mir <- mir
  df
}) %>%
  do.call(rbind, .)

n_sig <- ifelse(is.null(go) || nrow(go) == 0, 0, nrow(go))
cat("Significant terms at FDR<10%:", n_sig, "\n")

# --- Generate figure ---
if (n_sig > 0) {
  # Network plot (original Rmd approach)
  cat("Generating network plot...\n")
  df <- subset(go, mir == "hsa-mir-4532")[, c("Term", "Genes", "db", "mir")]
  df$db <- factor(df$db, c("GO_Molecular_Function_2025",
                            "GO_Biological_Process_2025",
                            "WikiPathways_2024_Human"))
  df$mir <- "hsa-miR-4532"

  p1 <- plot_mirna_network(df, shorten_factor = 1.5,
                            gene_node_size = 17,
                            mirna_node_size = 30,
                            legend_node_size = 6) +
    theme(plot.margin = margin(t = 10, r = 40, b = 10, l = 10),
          legend.box.margin = margin(10, 10, 10, 10),
          legend.margin = margin(20, 20, 20, 20)) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15))

  # Write PDF to a tempfile first (works everywhere: some filesystems, e.g.
  # S3 mounts, don't support PDF device writes), then copy to results dir
  out_pdf <- here("results", "Fig2D_enrichr_network.pdf")
  tmp_pdf <- tempfile(fileext = ".pdf")
  ggsave(tmp_pdf, p1, width = 10, height = 4)
  system2("cp", args = c(tmp_pdf, out_pdf))
  cat("Figure 2D saved to results/Fig2D_enrichr_network.pdf\n")

} else {
  # Fallback: dotplot of top nominal terms
  cat("No terms survive FDR<10%. Generating dotplot fallback...\n")

  # Collect top nominal terms from each database
  top_terms <- lapply(names(result), function(i) {
    res <- result[[i]]
    if (nrow(res) == 0) return(NULL)
    res %>%
      mutate(db = i) %>%
      arrange(P.value) %>%
      head(10)
  }) %>%
    do.call(rbind, .)

  if (is.null(top_terms) || nrow(top_terms) == 0) {
    cat("No enrichR results at all. Skipping Figure 6.\n")
    quit(save = "no", status = 0)
  }

  # Clean term names for display
  top_terms$Term_short <- sapply(strsplit(top_terms$Term, " \\("), function(i) i[[1]])

  p <- top_terms %>%
    mutate(Term_short = fct_reorder(Term_short, -log10(P.value))) %>%
    ggplot(aes(x = db, y = Term_short, size = -log10(P.value), color = -log10(P.value))) +
    geom_point() +
    scale_color_gradient(low = "steelblue", high = "firebrick",
                         name = expression(-log[10](p))) +
    scale_size_continuous(name = expression(-log[10](p)), range = c(2, 8)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          strip.text = element_text(size = 10)) +
    xlab("") +
    ylab("") +
    ggtitle("miR-4532 target enrichment (nominal p-values)")

  ggsave(here("results", "Fig2D_enrichr_dotplot.png"), p, height = 6, width = 8, dpi = 300)
  cat("Figure 2D (dotplot fallback) saved to results/Fig2D_enrichr_dotplot.png\n")
}
