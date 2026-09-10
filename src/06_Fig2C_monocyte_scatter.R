#!/usr/bin/env Rscript
# 06_Fig2C_monocyte_scatter.R
# Figure 2C: Monocyte × miR-4532 scatterplot per group
# Pearson correlation, ggscatter, faceted by group with n labels
# Output: results/Fig2C_monocyte_scatter.pdf
#
# Usage:
#   Rscript src/06_Fig2C_monocyte_scatter.R

library(here)
library(tidyverse)
library(ggpubr)

# --- Load data ---
cat("Loading data...\n")
mirna    <- readRDS(here("results", "01_mirna_filtered.rds"))
clinical <- readRDS(here("results", "01_clinical.rds"))

# --- Prepare data ---
cbc_df <- cbind(clinical[, c("Monocyte", "Group")],
               mir4532 = as.numeric(mirna["20518933_hsa-miR-4532", ])) %>%
  gather(cc, value, -c("Group", "mir4532")) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(
    value = as.numeric(value),
    Group = str_to_lower(Group)
  )

# Add group sample sizes
cbc_df$Group[cbc_df$Group %in% "control"] <- "control (n=25)"
cbc_df$Group[cbc_df$Group %in% "ambulatory"] <- "ambulatory (n=20)"
cbc_df$Group[cbc_df$Group %in% "non-ambulatory"] <- "non-ambulatory (n=16)"
cbc_df <- cbc_df %>%
  mutate(Group = factor(Group, c("control (n=25)", "ambulatory (n=20)",
                                  "non-ambulatory (n=16)")))

# --- Generate scatterplot ---
cat("Generating scatterplot...\n")
p1 <- cbc_df %>%
  ggscatter(x = "value", y = "mir4532",
            color = "Group",
            palette = "jco",
            add = "reg.line",
            conf.int = TRUE,
            cor.coef = TRUE,
            cor.method = "spearman",
            ylab = "miR-4532 expression",
            xlab = "Monocyte proportion",
            facet.by = "Group",
            scales = "free")

# Write PDF to a tempfile first (works everywhere: some filesystems, e.g. S3
# mounts, don't support PDF device writes), then copy to results directory
out_pdf <- here("results", "Fig2C_monocyte_scatter.pdf")
tmp_pdf <- tempfile(fileext = ".pdf")
ggsave(tmp_pdf, p1, height = 3, width = 6)
system2("cp", args = c(tmp_pdf, out_pdf))

cat("Figure 2C saved to results/Fig2C_monocyte_scatter.pdf\n")
