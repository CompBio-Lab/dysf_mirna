#!/usr/bin/env Rscript
# 04_Fig2A_mir4532boxplot.R
# Figure2A: miR-4532 boxplot by group
# Two-stage candidate selection: NonAmb vs Ctrl FDR<10%, re-test in Amb vs NonAmb
# Output: results/Fig2A_mir4532boxplot.png
#
# Usage:
#   Rscript src/04_Figure2A.R

library(here)
library(tidyverse)
library(ggpubr)

# --- Load data ---
cat("Loading data...\n")
mirna    <- readRDS(here("results", "01_mirna_filtered.rds"))
clinical <- readRDS(here("results", "01_clinical.rds"))
top1     <- readRDS(here("results", "03_de_results.rds"))

# --- Two-stage candidate selection (exact Rmd logic) ---
# Stage 1: Select miRNAs with FDR < 10% in control vs. non-ambulatory
sel_mirs <- top1 %>%
  filter(contrast == "control vs. non-ambulatory") %>%
  filter(adj.P.Val < 0.1) %>%
  pull(feature)

cat("Candidate miRNAs from control vs. non-ambulatory (FDR<10%):", length(sel_mirs), "\n")
cat("  ", sel_mirs, "\n")

# Stage 2: Re-test in ambulatory vs. non-ambulatory with restricted BH correction
sig <- top1 %>%
  filter(feature %in% sel_mirs) %>%
  filter(contrast == "ambulatory vs. non-ambulatory") %>%
  mutate(adj.P.Val = p.adjust(P.Value, "BH")) %>%
  filter(adj.P.Val < 0.1)

cat("Confirmed miRNAs (restricted FDR<10% in amb vs. non-amb):", nrow(sig), "\n")
if (nrow(sig) > 0) print(sig[, c("feature", "logFC", "P.Value", "adj.P.Val")])

# --- Boxplot for miR-4532 ---
cat("Generating miR-4532 boxplot...\n")

df <- t(mirna[sapply(c("hsa-miR-4532"), function(i) grep(i, rownames(mirna))), ]) %>%
  as.data.frame() %>%
  mutate(group = clinical$Group) %>%
  gather(mirna, exp, -group) %>%
  mutate(mirna = sapply(strsplit(mirna, "_"), function(i) i[[2]]))

df$mirna <- gsub("hsa-", "", df$mirna)
df$group <- as.character(df$group)
df$group[df$group %in% "Ambulatory"] <- "ambulatory"
df$group[df$group %in% "Control"] <- "control"
df$group[df$group %in% "Non-ambulatory"] <- "non-ambulatory"
df$group <- factor(df$group, c("control", "ambulatory", "non-ambulatory"))

p <- ggboxplot(df, x = "group", y = "exp", color = "group",
               add = "jitter",
               palette = c("#999999", "#E69F00", "#56B4E9"))
p2 <- facet(p, "mirna", scales = "free") +
  geom_pwc(method = "t.test") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ylab(expression("log"[2]~"expression")) +
  theme_classic() +
  theme(strip.text = element_text(size = 25),
        legend.position = "none")

ggsave(here("results", "Fig2A_mir4532boxplot.png"), p2, height = 3, width = 3, dpi = 300)
cat("Figure 2A saved to results/Fig2A_mir4532boxplot.png\n")
