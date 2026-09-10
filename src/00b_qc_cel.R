#!/usr/bin/env Rscript
# 00b_qc_cel.R
# Quality control report for raw CEL files (Affymetrix miRNA 4.0).
#
# QC levels:
#   1. Raw intensity distributions (density + boxplot)
#   2. RLE / NUSE from probe-level model (classic Affymetrix array metrics)
#   3. Spike-in control ladder (miRNA 4.0 "spike_in-control-N" probes:
#      per-array rank correlation vs the consensus ladder + mean-signal
#      outlier check; self-calibrating, no concentration map needed)
#   4. Post-RMA distributions, PCA, and sample-correlation heatmap
#
# Inputs:  cel_files/*.CEL
# Outputs: results/qc/cel_qc_report.pdf   (multi-page QC report)
#          results/qc/qc_*.png           (each plot as an individual PNG)
#          results/qc/cel_qc_metrics.csv (per-array metrics + PASS/REVIEW flag)
#          results/qc/cel_qc_summary.csv (metric ranges, thresholds, failures)
#
# Flagging heuristics (REVIEW if any fail):
#   NUSE median > 1.05 | RLE |median| > 0.2 | RLE IQR > 0.5 | spike-in check fail
# Exclude an array only if flagged by >= 2 independent metrics.
#
# Uses only oligo + pd.mirna.4.0 (already installed by the Dockerfile).
#
# Usage:
#   Rscript src/00b_qc_cel.R

library(here)
library(oligo)
library(pd.mirna.4.0)

NUSE_MED_MAX <- 1.05
RLE_MED_MAX  <- 0.2
RLE_IQR_MAX  <- 0.5

out_dir <- here("results", "qc")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load CEL files ---
cel_files <- list.files(here("cel_files"), pattern = "\\.CEL$",
                        full.names = TRUE, ignore.case = TRUE)
if (length(cel_files) == 0)
  stop("No CEL files found in cel_files/. Place them there first.")
cat("Found", length(cel_files), "CEL files\n")

raw <- read.celfiles(cel_files, pkgname = "pd.mirna.4.0", checkType = FALSE)
# Sample labels: 4th underscore-delimited field (project convention)
sn <- sapply(strsplit(sampleNames(raw), "_"), function(i) i[[4]])
sn[is.na(sn)] <- sampleNames(raw)[is.na(sn)]
sampleNames(raw) <- make.unique(sn)

# --- RMA (needed for post-normalization QC + spike-in check) ---
cat("Running RMA for post-normalization QC...\n")
eset <- rma(raw)
expr_mat <- exprs(eset)

# --- RLE / NUSE via probe-level model ---
cat("Fitting probe-level model for RLE/NUSE (this takes a few minutes)...\n")
plm <- fitProbeLevelModel(raw)
rle_vals <- RLE(plm, type = "values")    # probesets x arrays, centered at 0
nuse_vals <- NUSE(plm, type = "values")  # probesets x arrays, centered at ~1
colnames(rle_vals) <- colnames(nuse_vals) <- sampleNames(raw)

rle_med  <- apply(rle_vals, 2, median)
rle_iqr  <- apply(rle_vals, 2, IQR)
nuse_med <- apply(nuse_vals, 2, median)

# --- Spike-in control ladder check ---
# miRNA 4.0 uses "spike_in-control-N" probes (not the classic poly-A
# bioB/bioC/bioD/creX set). Without a concentration map, use a
# self-calibrating check: detected spike-ins (cross-array mean > 1) form a
# consensus ladder; each array must (a) reproduce the ladder rank order
# (Spearman rho >= 0.9) and (b) have a mean spike-in signal within 2 SD of
# the cross-array mean.
cat("Checking spike-in controls...\n")
ann <- read.csv(here("data", "mirna", "miRNA-4_0-st-v1.annotations.20160922.csv"),
                row.names = 1)
spike <- ann[grepl("spike", ann$Sequence.Type, ignore.case = TRUE), ]
spike_present <- rownames(expr_mat) %in% spike$Probe.Set.Name
spike_mat <- expr_mat[spike_present, , drop = FALSE]
spike_ok <- rep(NA, ncol(expr_mat))
spike_rho <- rep(NA_real_, ncol(expr_mat))
detected <- rowMeans(spike_mat) > 1
if (sum(detected) >= 3) {
  spike_det <- spike_mat[detected, , drop = FALSE]
  consensus <- rowMeans(spike_det)
  cat("Detected spike-in ladder (", nrow(spike_det), " probes): ",
      paste(round(sort(consensus), 2), collapse = " < "), "\n", sep = "")
  spike_rho <- apply(spike_det, 2, function(x) cor(x, consensus, method = "spearman"))
  spike_mean <- colMeans(spike_det)
  spike_out <- abs(spike_mean - mean(spike_mean)) > 2 * sd(spike_mean)
  spike_ok <- (spike_rho >= 0.9) & !spike_out
} else {
  cat("WARNING: fewer than 3 detected spike-in controls; skipping ladder check.\n")
}

# --- Per-array metrics + flags ---
metrics <- data.frame(
  sample = sampleNames(raw),
  # NB: oligo:: qualifier required — lubridate (tidyverse) masks pm()
  mean_raw_log2 = round(colMeans(log2(oligo::pm(raw))), 3),
  rle_median = round(rle_med, 4),
  rle_iqr = round(rle_iqr, 4),
  nuse_median = round(nuse_med, 4),
  spikein_ladder_rho = round(spike_rho, 4),
  spikein_order_pass = spike_ok,
  stringsAsFactors = FALSE
)
metrics$flag_nuse <- metrics$nuse_median > NUSE_MED_MAX
metrics$flag_rle  <- abs(metrics$rle_median) > RLE_MED_MAX | metrics$rle_iqr > RLE_IQR_MAX
metrics$flag_spike <- !is.na(metrics$spikein_order_pass) & !metrics$spikein_order_pass
metrics$n_flags <- rowSums(metrics[, c("flag_nuse", "flag_rle", "flag_spike")])
metrics$qc_call <- ifelse(metrics$n_flags == 0, "PASS", "REVIEW")
flagged <- metrics$sample[metrics$qc_call == "REVIEW"]
flag_col <- ifelse(metrics$qc_call == "REVIEW", "#D55E00", "grey30")

write.csv(metrics, file.path(out_dir, "cel_qc_metrics.csv"), row.names = FALSE)

# --- Summary table (metric ranges, thresholds, failure counts) ---
fmt_range <- function(x, digits = 3)
  paste(round(min(x, na.rm = TRUE), digits), "to", round(max(x, na.rm = TRUE), digits))
rle_med_fail  <- metrics$sample[abs(metrics$rle_median) > RLE_MED_MAX]
rle_iqr_fail  <- metrics$sample[metrics$rle_iqr > RLE_IQR_MAX]
nuse_fail     <- metrics$sample[metrics$nuse_median > NUSE_MED_MAX]
rho_fail      <- metrics$sample[!is.na(metrics$spikein_ladder_rho) &
                                metrics$spikein_ladder_rho < 0.9]
spk_mean_fail <- metrics$sample[metrics$flag_spike & !metrics$sample %in% rho_fail]
summary_tab <- data.frame(
  metric = c("RLE median", "RLE IQR", "NUSE median", "Spike-in ladder rho",
             "Spike-in mean signal", "Overall PASS", "Overall REVIEW"),
  observed_range = c(fmt_range(metrics$rle_median), fmt_range(metrics$rle_iqr),
                     fmt_range(metrics$nuse_median), fmt_range(metrics$spikein_ladder_rho),
                     "within 2 SD of cross-array mean", "", ""),
  threshold = c("|x| > 0.2", "> 0.5", "> 1.05", "< 0.9",
                "> 2 SD from cross-array mean", "", ""),
  n_flagged = c(length(rle_med_fail), length(rle_iqr_fail), length(nuse_fail),
                length(rho_fail), length(spk_mean_fail),
                sum(metrics$qc_call == "PASS"), sum(metrics$qc_call == "REVIEW")),
  samples = c(paste(rle_med_fail, collapse = "; "), paste(rle_iqr_fail, collapse = "; "),
              paste(nuse_fail, collapse = "; "), paste(rho_fail, collapse = "; "),
              paste(spk_mean_fail, collapse = "; "), "",
              paste(metrics$sample[metrics$qc_call == "REVIEW"], collapse = "; ")),
  stringsAsFactors = FALSE
)
write.csv(summary_tab, file.path(out_dir, "cel_qc_summary.csv"), row.names = FALSE)

# --- PCA (used by the PCA plot; computed once) ---
pca <- prcomp(t(expr_mat), center = TRUE, scale. = TRUE, rank. = 5)
ve <- round(summary(pca)$importance[2, 1:5] * 100, 1)

# --- Plot definitions ---
# Each plot is a function so it can be rendered identically twice:
# once as an individual PNG and once as a page of the combined PDF.
plots <- list(
  raw_density = function() {
    par(mar = c(5, 4.5, 3, 1))
    oligo::hist(raw, main = "Raw intensity distributions",
                xlab = "log2 intensity", col = flag_col, lty = 1)
    legend("topright", legend = c("PASS", "REVIEW"),
           col = c("grey30", "#D55E00"), lty = 1, cex = 0.8, bty = "n")
  },
  raw_boxplot = function() {
    par(mar = c(9, 4.5, 3, 1))
    oligo::boxplot(raw, main = "Raw intensity boxplots", las = 2,
                   cex.axis = 0.5, col = flag_col, ylab = "log2 intensity")
  },
  rle = function() {
    par(mar = c(9, 4.5, 3, 1))
    boxplot(rle_vals, las = 2, cex.axis = 0.5, col = flag_col, outline = FALSE,
            main = "RLE (Relative Log Expression)", ylab = "RLE")
    abline(h = c(-RLE_IQR_MAX / 2, RLE_IQR_MAX / 2), lty = 2, col = "grey50")
  },
  nuse = function() {
    par(mar = c(9, 4.5, 3, 1))
    boxplot(nuse_vals, las = 2, cex.axis = 0.5, col = flag_col, outline = FALSE,
            main = "NUSE (Normalized Unscaled Standard Error)", ylab = "NUSE")
    abline(h = NUSE_MED_MAX, lty = 2, col = "#D55E00")
  },
  rma_boxplot = function() {
    par(mar = c(9, 4.5, 3, 1))
    boxplot(expr_mat, las = 2, cex.axis = 0.5, col = flag_col, outline = FALSE,
            main = "RMA-normalized expression", ylab = "log2 expression")
  },
  rma_density = function() {
    par(mar = c(5, 4.5, 3, 1))
    plot(density(expr_mat[, 1]), col = flag_col[1],
         main = "RMA-normalized distributions", xlab = "log2 expression",
         ylim = c(0, max(sapply(1:ncol(expr_mat),
                                function(i) max(density(expr_mat[, i])$y)))))
    for (i in 2:ncol(expr_mat))
      lines(density(expr_mat[, i]), col = flag_col[i])
  },
  pca = function() {
    par(mar = c(5, 4.5, 3, 1))
    plot(pca$x[, 1], pca$x[, 2], pch = 19, col = flag_col,
         xlab = sprintf("PC1 (%.1f%%)", ve[1]),
         ylab = sprintf("PC2 (%.1f%%)", ve[2]),
         main = "PCA of RMA-normalized data (all probe sets)")
    text(pca$x[, 1], pca$x[, 2], labels = sampleNames(raw), pos = 3, cex = 0.45)
  },
  cor_heatmap = function() {
    par(mar = c(7, 7, 3, 1))
    heatmap(cor(expr_mat), symm = TRUE, margins = c(6, 6),
            main = "Sample-sample Pearson correlation",
            col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(50))
  }
)
plot_array_image <- function(s)
  oligo::image(raw[, which(sampleNames(raw) == s)], main = paste("Array image:", s))

# --- Save each plot as an individual PNG in results/qc/ ---
cat("Saving individual QC plots to", out_dir, "...\n")
for (nm in names(plots)) {
  png(file.path(out_dir, paste0("qc_", nm, ".png")),
      width = 11, height = 8.5, units = "in", res = 300)
  plots[[nm]]()
  dev.off()
}
for (s in flagged) {
  png(file.path(out_dir, paste0("qc_array_image_", s, ".png")),
      width = 8, height = 8, units = "in", res = 300)
  plot_array_image(s)
  dev.off()
}

# --- Combined multi-page PDF (tempfile + cp: pdf() cannot write S3 mounts) ---
tmp_pdf <- tempfile(fileext = ".pdf")
pdf(tmp_pdf, width = 11, height = 8.5)
for (p in plots) p()
for (s in flagged) plot_array_image(s)
dev.off()

cp_ok <- system2("cp", c(tmp_pdf, file.path(out_dir, "cel_qc_report.pdf")))
if (cp_ok != 0) stop("Failed to copy QC PDF to ", out_dir)

# --- Summary ---
cat("\n======== QC summary ========\n")
cat("Arrays:", nrow(metrics), "| PASS:", sum(metrics$qc_call == "PASS"),
    "| REVIEW:", sum(metrics$qc_call == "REVIEW"), "\n")
if (length(flagged) > 0) {
  cat("Flagged arrays:\n")
  print(metrics[metrics$qc_call == "REVIEW",
                c("sample", "rle_median", "rle_iqr", "nuse_median",
                  "spikein_order_pass", "n_flags")])
}
cat("Report:  ", file.path(out_dir, "cel_qc_report.pdf"), "\n")
cat("Plots:   ", length(plots) + length(flagged), "PNGs in", out_dir, "\n")
cat("Metrics: ", file.path(out_dir, "cel_qc_metrics.csv"), "\n")
cat("Summary: ", file.path(out_dir, "cel_qc_summary.csv"), "\n")
