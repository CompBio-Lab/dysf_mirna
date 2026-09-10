#!/usr/bin/env Rscript
# 00_cel_to_normalized.R
# Preprocessing: CEL files → RMA-normalized expression matrix (oligo)
# Output: data/mirna/normalized_data_celderived.TXT  (PRIMARY analysis matrix)
# NOTE: the legacy Affymetrix TAC export (normalized_data.TXT) is kept as a
# separate file and is never overwritten by this script.
#
# Prerequisites:
#   - Place 74 CEL files in the cel_files/ directory (mounted separately)
#   - The pd.mirna.4.0 package must be installed (included in Docker image)
#
# Output format matches the original normalized_data.TXT (Affymetrix TAC export):
#   - Row names: numeric probe set IDs (e.g. "20518933"). oligo/rma returns
#     manufacturer probe set names (e.g. "MIMAT0019071_st"); these are mapped
#     to numeric IDs via the pd.mirna.4.0 featureSet table (fsetid column),
#     which is identical to TAC's "Probe Set ID" numbering.
#   - Column names: full CEL filenames (e.g. "060815_Proof_MD_MD00101P_miRNA_4_1_A01.CEL")
#     — script 01 extracts the sample ID as the 4th underscore-delimited field.
#     (TAC appends ".rma-dabg-Signal" instead of ".CEL"; both sit on the last
#     field, so downstream extraction is unaffected.)
#   - All 36,353 probe sets on the miRNA-4.0 array are written (human, other
#     species, controls). Script 01 filters to human expressed probes.
#
# Usage:
#   Rscript src/00_cel_to_normalized.R

library(here)
library(oligo)
library(DBI)

# --- Configuration ---
cel_dir   <- here("cel_files")
out_file  <- here("data", "mirna", "normalized_data_celderived.TXT")

# --- Read CEL files ---
cat("Reading CEL files from:", cel_dir, "\n")
cel_files <- list.celfiles(cel_dir, full.names = TRUE)
cat("Found", length(cel_files), "CEL files\n")

if (length(cel_files) == 0) {
  stop("No CEL files found in cel_files/. Please mount the CEL directory.")
}

# Read with pd.mirna.4.0 platform design package
# checkType=FALSE (an argument of read.celfiles, NOT rma) allows miRNA-4_1
# arrays to be read with the 4.0 annotation package
cat("Reading CEL files with oligo...\n")
raw <- read.celfiles(cel_files, pkgname = "pd.mirna.4.0", checkType = FALSE)

# --- RMA normalization ---
cat("Running RMA normalization...\n")
eset <- rma(raw)
expr_mat <- exprs(eset)

# --- Map probe set names to numeric probe set IDs (TAC format) ---
# oligo returns man_fsetid (e.g. "MIMAT0019071_st"); TAC uses fsetid (20518933)
cat("Mapping probe set names to numeric probe set IDs...\n")
fmap <- dbGetQuery(db(pd.mirna.4.0), "SELECT fsetid, man_fsetid FROM featureSet")
idx <- match(rownames(expr_mat), fmap$man_fsetid)
if (any(is.na(idx))) {
  stop("Failed to map ", sum(is.na(idx)), " probe set names to numeric IDs")
}
rownames(expr_mat) <- fmap$fsetid[idx]

# --- Write output ---
# Keep full CEL filenames as column names; 01_data_import.R extracts sample IDs
cat("Writing normalized data to:", out_file, "\n")
cat("Dimensions:", nrow(expr_mat), "probes x", ncol(expr_mat), "samples\n")
write.table(expr_mat, file = out_file, sep = "\t", quote = FALSE,
            row.names = TRUE, col.names = NA)

cat("Done. Normalized data saved.\n")
