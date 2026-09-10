#!/usr/bin/env Rscript
# 05_Fig2B_cbc_heatmap.R
# Figure 2B: CBC × miR-4532 correlation heatmap per group
# Per-group Spearman correlations, ComplexHeatmap with |rho|>0.4 labels
# Output: results/Fig2B_cbc_heatmap.pdf
#
# Usage:
#   Rscript src/05_figure4_cbc_heatmap.R

library(here)
library(ComplexHeatmap)
library(circlize)
library(grid)

# --- Load data ---
cat("Loading data...\n")
mirna    <- readRDS(here("results", "01_mirna_filtered.rds"))
clinical <- readRDS(here("results", "01_clinical.rds"))
top_mirna_subset <- readRDS(here("results", "03_de_sig_features.rds"))

# --- Compute per-group correlations ---
cat("Computing per-group Spearman correlations...\n")

corlist <- lapply(levels(clinical$Group), function(i) {
  # Get significant miRNA probes
  mirna_top <- mirna[unique(as.character(sapply(top_mirna_subset$feature, function(j) {
    grep(j, rownames(mirna), value = TRUE)
  }))), clinical$Group == i]

  # Get CBC values
  cbc_top <- apply(clinical[clinical$Group == i, c("WBC", "RBC", "Hemoglobin",
                     "Hematocrit", "MCV", "RDW", "Platelet", "MPV",
                     "Neut", "Lymph", "Monocyte", "Eosinophils", "Basophils",
                     "Immature.Granulocytes")], 2, as.numeric)
  rownames(cbc_top) <- rownames(clinical[clinical$Group == i, ])
  cbc_top <- t(cbc_top)

  # Correlate CBC with miR-4532 only
  cormat <- cor(t(cbc_top), t(mirna_top[c("20518933_hsa-miR-4532"), , drop = FALSE]),
                use = "pairwise.complete.obs")
  cormat
})

cordat <- do.call(cbind, corlist)
colnames(cordat) <- NULL

# --- Column annotation ---
group <- c("control", "ambulatory", "non-ambulatory")
mir <- rep(c("miR-4532"), 3)
col_ha <- columnAnnotation(
  Group = factor(group, group),
  col = list(Group = c("control" = "#999999", "ambulatory" = "#E69F00",
                        "non-ambulatory" = "#56B4E9"))
)

# --- Draw heatmap ---
# Write PDF to a tempfile first (works everywhere: some filesystems, e.g. S3
# mounts, don't support pdf() device writes), then copy to results directory
cat("Drawing heatmap...\n")
out_pdf <- here("results", "Fig2B_cbc_heatmap.pdf")
tmp_pdf <- tempfile(fileext = ".pdf")
pdf(tmp_pdf, height = 3)
Heatmap(cordat, cluster_columns = FALSE,
        top_annotation = col_ha,
        column_split = mir,
        name = "Correlation",
        border = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          val <- cordat[i, j]
          if (abs(val) > 0.4) {
            grid.text(sprintf("%.2f", val), x, y, gp = gpar(fontsize = 15))
          }
        })
dev.off()
system2("cp", args = c(tmp_pdf, out_pdf))

cat("Figure 2B saved to results/Fig2B_cbc_heatmap.pdf\n")
