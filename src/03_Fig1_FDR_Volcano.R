#!/usr/bin/env Rscript
# 03_Fig1_fdr_volcano.R
# Figure 1: FDR rank plot (A) + volcano plots (B) for 3 contrasts
# Uses limma with ~Group + sex + age + 0 (no batch covariate)
# Output: results/Fig1_fdr_volcano.png
#
# Usage:
#   Rscript src/03_figure2_volcano.R

library(here)
library(tidyverse)
library(limma)
library(cowplot)
library(ggrepel)

# --- Load filtered data ---
cat("Loading filtered data...\n")
mirna    <- readRDS(here("results", "01_mirna_filtered.rds"))
clinical <- readRDS(here("results", "01_clinical.rds"))

# --- Design matrix (no batch covariate, matching Rmd) ---
design <- model.matrix(~Group + sex + age + 0, data = clinical)
colnames(design) <- gsub("Group", "", colnames(design))
colnames(design) <- gsub("n-a", "nA", colnames(design))

cat("Design matrix columns:", colnames(design), "\n")

# --- Fit limma ---
lmfit <- lmFit(mirna, design)
cont <- makeContrasts(Ambulatory - Control, NonAmbulatory - Ambulatory,
                      NonAmbulatory - Control, levels = design)
lmfit.cont <- contrasts.fit(lmfit, cont)
lmfit.cont.ebayes <- eBayes(lmfit.cont)

# --- Extract top tables ---
top1 <- lapply(colnames(cont), function(contrast) {
  topTable(lmfit.cont.ebayes, coef = contrast,
           adjust.method = "BH", n = nrow(lmfit.cont.ebayes)) %>%
    as.data.frame() %>%
    mutate(contrast = contrast) %>%
    mutate(feature = rownames(.))
}) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(feature = sapply(strsplit(feature, "_"), function(i) i[[2]])) %>%
  group_by(contrast) %>%
  mutate(n = 1:n())

# --- Rename contrasts ---
top1$contrast[top1$contrast == "Ambulatory - Control"] <- "control vs. ambulatory"
top1$contrast[top1$contrast == "NonAmbulatory - Control"] <- "control vs. non-ambulatory"
top1$contrast[top1$contrast == "NonAmbulatory - Ambulatory"] <- "ambulatory vs. non-ambulatory"
top1$contrast <- factor(top1$contrast,
                         c("control vs. ambulatory", "control vs. non-ambulatory",
                           "ambulatory vs. non-ambulatory"))

# --- Panel A: FDR rank plot ---
p1 <- top1 %>%
  mutate(logp = -log10(P.Value),
         FC = ifelse(logFC > 0, "UP", "DOWN")) %>%
  ggplot(aes(x = n, y = adj.P.Val, color = contrast)) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  theme_bw() +
  xlab("Number of significant miRNAs") +
  ylab("BH-FDR") +
  theme(legend.position = "top") +
  guides(color = guide_legend(ncol = 1)) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  annotate("text", x = 130, y = 0.15, label = "BH-FDR=10%", size = 3)

# --- Panel B: Volcano plots ---
top_mirna <- top1 %>%
  mutate(logp = -log10(P.Value),
         sig = ifelse(adj.P.Val < 0.1, "Significant", "Not-significant")) %>%
  mutate(feature = gsub("hsa-", "", feature))

top_mirna_subset <- subset(top_mirna, adj.P.Val < 0.1)
top_mirna_subset2 <- top_mirna %>%
  dplyr::group_by(contrast) %>%
  dplyr::slice(1) %>%
  dplyr::filter(adj.P.Val > 0.1)

## save mirna (non-ambulatory vs. control)
top_mirna_subset %>% 
  filter(contrast == "control vs. non-ambulatory") %>% 
  dplyr::select(feature, logFC, AveExpr, P.Value, adj.P.Val) %>% 
  mutate(logFC = signif(logFC, 2),
         AveExpr = signif(AveExpr, 2),
         P.Value = signif(P.Value, 2),
         adj.P.Val = signif(adj.P.Val, 2)) %>% 
  write.csv(here("results", "amb_vs_control_sig_mirs.csv"))

p2 <- top_mirna %>%
  ggplot(aes(x = logFC, y = logp, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point() +
  facet_wrap(~contrast, scales = "free") +
  ggrepel::geom_text_repel(data = top_mirna_subset2,
                           aes(label = feature), color = "black",
                           size = 5, max.overlaps = 50) +
  ggrepel::geom_text_repel(data = top_mirna_subset,
                           aes(label = feature), size = 5, max.overlaps = 50) +
  scale_color_manual(
    values = c("Not-significant" = "darkgrey",
               "Significant"     = "blue"),
    breaks = c("Significant", "Not-significant")) +
  theme_classic() +
  theme(legend.position = c(0.08, 0.85),
        strip.text = element_text(size = 16),
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  xlab(expression("log"[2]~"fold-change")) +
  ylab(expression("-log"[10]~"P-value"))

# --- Combine panels ---
p3 <- cowplot::plot_grid(p1, p2, labels = c("A", "B"), rel_widths = c(1, 2.2))

ggsave(here("results", "Fig1_fdr_volcano.png"), p3, height = 5, width = 15, dpi = 300)
cat("Figure 1 saved to results/Fig1_fdr_volcano.png\n")

# --- Save DE results for downstream scripts ---
saveRDS(top1, here("results", "03_de_results.rds"))
saveRDS(top_mirna_subset, here("results", "03_de_sig_features.rds"))
cat("DE results saved for downstream use.\n")
