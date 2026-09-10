#!/usr/bin/env Rscript
# 02_FigS1_pca_heatmap.R
# Figure S1: PCA × clinical variable heatmap
# Uses pcaHeatmap() from utils.R with ANOVA p-values
# Output: results/FigS1_pca_heatmap.png
#
# Usage:
#   Rscript src/02_figure1_pca_heatmap.R

library(here)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)

source(here("src", "utils.R"))

# --- Load filtered data ---
cat("Loading filtered data...\n")
mirna    <- readRDS(here("results", "01_mirna_filtered.rds"))
clinical <- readRDS(here("results", "01_clinical.rds"))

# --- PCA (rank=5, matching Rmd) ---
cat("Running PCA (rank=5)...\n")
pcs <- prcomp(t(mirna), scale. = TRUE, center = TRUE, rank. = 5)
cat("Variance explained:\n")
print(summary(pcs))

# Label PCs with % variance explained
# (sdev may include all components even with rank.=5 — use full sdev for total, keep first k)
var_pct <- round(100 * pcs$sdev^2 / sum(pcs$sdev^2), 1)[seq_len(ncol(pcs$x))]
colnames(pcs$x) <- paste0("PC", seq_len(ncol(pcs$x)), " (", var_pct, "%)")

# --- Prepare clinical variables for pcaHeatmap ---
clinical$Sex <- clinical$sex
clinical$Age <- clinical$age

# --- Generate heatmap ---
cat("Generating PCA heatmap...\n")
p0 <- pcaHeatmap(pcs$x, demo = clinical[, c("Group", "RNA.Extraction.Batch", "Sex", "Age")])

# --- Save ---
ggsave(here("results", "FigS1_pca_heatmap.png"), p0, height = 3, width = 5, dpi = 300)
cat("FigS1_pca_heatmap saved to results/FigS1_pca_heatmap.png\n")
