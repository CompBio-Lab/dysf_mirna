library(tidygraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(grid) # unit()

pcaHeatmap = function(pcs, demo) {
  pvalheatmap <- matrix(0, ncol = ncol(demo), nrow = ncol(pcs))
  rownames(pvalheatmap) <- colnames(pcs)
  colnames(pvalheatmap) <- colnames(demo)
  for (i in 1:ncol(pcs)) {
    for (j in 1:ncol(demo)) {
      pvalheatmap[i, j] <- summary(aov(lm(pcs[, i] ~ demo[, 
                                                          j])))[[1]][1, "Pr(>F)"]
    }
  }
  pvalheatmap[pvalheatmap < 0.01] <- 0.01
  pvalheatmap[pvalheatmap > 0.1] <- 1
  pvalheatmap[pvalheatmap > 0.01 & pvalheatmap < 0.05] <- 0.05
  pvalheatmap[pvalheatmap > 0.05 & pvalheatmap < 0.1] <- 0.1
  pvalheatmap[pvalheatmap == "0.01"] <- "p < 0.01"
  pvalheatmap[pvalheatmap == "0.05"] <- "0.01 < p < 0.05"
  pvalheatmap[pvalheatmap == "0.1"] <- "0.05 < p < 0.10"
  pvalheatmap[pvalheatmap == "1"] <- "p > 0.10"
  pvalheatmap %>% 
    as.data.frame %>% 
    mutate(Variable = rownames(.)) %>% 
    tidyr::gather(Threshold, Value, -Variable) %>% 
    mutate(Threshold = factor(Threshold, levels = unique(Threshold))) %>% 
    mutate(Value = factor(Value, levels = c("p < 0.01", "0.01 < p < 0.05", "0.05 < p < 0.10", "p > 0.10"))) %>% 
    ggplot(aes(Threshold, Variable)) + 
    geom_tile(aes(fill = Value), colour = "white") + 
    scale_fill_manual(values = rev(RColorBrewer::brewer.pal(n = 8, name = "Blues")[c(2, 4, 6, 8)])) + 
    customTheme(sizeStripFont = 10, xAngle = 40, hjust = 1, vjust = 1, xSize = 10, ySize = 10, xAxisSize = 10, yAxisSize = 10) + 
    xlab("") + 
    ylab("") +
    labs(fill = "Significance intervals")
}


customTheme = function (sizeStripFont, xAngle, hjust, vjust, xSize, ySize, 
          xAxisSize, yAxisSize) {
  theme(strip.background = element_rect(colour = "black", fill = "white", 
                                        size = 1), strip.text.x = element_text(size = sizeStripFont), 
        strip.text.y = element_text(size = sizeStripFont), axis.text.x = element_text(angle = xAngle, 
                                                                                      hjust = hjust, vjust = vjust, size = xSize, color = "black"), 
        axis.text.y = element_text(size = ySize, color = "black"), 
        axis.title.x = element_text(size = xAxisSize, color = "black"), 
        axis.title.y = element_text(size = yAxisSize, color = "black"), 
        panel.background = element_rect(fill = "white", color = "black"))
}

plot_mirna_network <- function(
    df,
    gene_node_size = 6,
    mirna_node_size = 8,
    legend_node_size = 4,   # size of node type symbols in legend
    label_size = 3,
    arrow_mm = 3,
    edge_width = 0.6,
    edge_alpha = 0.7,
    shorten_factor = 1.8
) {
  # Expand genes
  edges_long <- df %>%
    mutate(Gene = str_split(Genes, ";")) %>%
    select(Term, mir, Gene) %>%
    unnest(Gene) %>%
    mutate(Gene = trimws(Gene))
  
  # Directed edge list: miRNA -> gene
  edges <- edges_long %>%
    transmute(from = mir, to = Gene, Term)
  
  # Node table
  mir_set <- unique(df$mir)
  nodes <- tibble(name = unique(c(edges$from, edges$to))) %>%
    mutate(
      type = if_else(name %in% mir_set, "miRNA", "gene"),
      Term = if_else(type == "gene",
                     edges_long$Term[match(name, edges_long$Gene)],
                     NA_character_)
    )
  
  # GO term factor
  term_levels <- nodes %>% filter(type == "gene") %>% distinct(Term) %>% pull(Term)
  nodes <- nodes %>%
    mutate(term_fill = factor(if_else(type == "gene", Term, NA_character_),
                              levels = term_levels))
  
  # Build graph
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  # Colors
  base_cols <- c("#6bb5ff", "#8ee091", "#f6b26b", "#d5a6bd", "#ffd966",
                 "#b6d7a8", "#a4c2f4", "#e06666", "#93c47d", "#8e7cc3")
  term_palette <- setNames(base_cols[seq_along(term_levels)], term_levels)
  
  # Caps
  start_cap_pts <- mirna_node_size * shorten_factor
  end_cap_pts   <- gene_node_size * shorten_factor
  
  # Plot
  set.seed(42)
  g <- ggraph(graph, layout = "fr") +
    geom_edge_link(
      arrow = arrow(length = unit(arrow_mm, "mm"), type = "closed"),
      start_cap = ggraph::circle(start_cap_pts, "pt"),
      end_cap   = ggraph::circle(end_cap_pts, "pt"),
      width = edge_width,
      alpha = edge_alpha,
      colour = "grey30"
    ) +
    geom_node_point(
      aes(fill = term_fill, shape = type,
          size = ifelse(type == "miRNA", mirna_node_size, gene_node_size)),
      colour = "black"
    ) +
    scale_size_identity() +
    geom_node_label(
      aes(label = name, fill = ifelse(type == "miRNA", "white", NA)),
      size = label_size,
      label.padding = unit(0.15, "lines"),
      label.size = NA,
      show.legend = FALSE
    ) +
    scale_shape_manual(values = c(miRNA = 23, gene = 21), name = "Node type") +
    scale_fill_manual(
      values = c(term_palette, "white" = "white"),
      breaks = term_levels,
      na.translate = FALSE,
      name = "GO term (genes)"
    ) +
    guides(
      fill  = guide_legend(override.aes = list(shape = 21, size = gene_node_size, colour = "black")),
      shape = guide_legend(override.aes = list(size = legend_node_size, fill = "grey90"))
    ) +
    theme_void() +
    ggtitle("miRNA -> Gene Network")
  
  return(g)
}
