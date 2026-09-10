#!/usr/bin/env Rscript
# 01_data_import.R
# Load normalized miRNA data + clinical/CBC data, filter probes, save RDS
# Output: data/mirna/mirna.csv, data/clinical/clinical.csv,
#         results/01_mirna_filtered.rds, results/01_clinical.rds
#
# This script reproduces the exact "import data" chunk from the Rmd pipeline.
# Primary input: normalized_data_celderived.TXT (oligo RMA from script 00).
# The legacy TAC export (normalized_data.TXT) can be used instead by setting
# input_file below — kept for sensitivity comparison only.
#
# Usage:
#   Rscript src/01_data_import.R

library(here)
library(tidyverse)

# --- Analysis scope ---
# TRUE  = mature miRNAs only (Sequence.Type == "miRNA"; excludes stem-loop
#         hairpins, snoRNAs, scaRNAs, rRNAs, spike-ins) -> 449 probes (TAC
#         input) / 438 probes (CEL-derived input from script 00)
# FALSE = all human-annotated probes (original Rmd behavior) -> 666 / 651
mature_only <- TRUE

# --- Load clinical data (TWO files joined) ---
cat("Loading clinical data...\n")
clinical0 <- read.csv(here("data", "clinical", "Samples sent to Scripps_2015.04.07_age_sex.csv"),
                       row.names = 3)
clinical0$Group[clinical0$Group == "Non-Ambulatory"] <- "Non-ambulatory"

cbc <- read.csv(here("data", "cbc", "md_cbc.csv"))
cbc$Study.ID <- gsub("MDN-", "MDN-0", cbc$Study.ID)  # CRITICAL: pad MDN IDs

clinical <- cbc %>%
  dplyr::rename(STUDY_ID = Study.ID) %>%
  inner_join(mutate(clinical0, id = rownames(clinical0)), by = "STUDY_ID")
rownames(clinical) <- clinical$id

# --- Load miRNA data ---
cat("Loading miRNA expression data...\n")
input_file <- "normalized_data_celderived.TXT"   # primary: oligo RMA (script 00)
# input_file <- "normalized_data.TXT"            # legacy: Affymetrix TAC export
mirna <- read.delim(here("data", "mirna", input_file), row.names = 1)
mirna_ann <- read.csv(here("data", "mirna", "miRNA-4_0-st-v1.annotations.20160922.csv"),
                       row.names = 1)
mirna_ann <- subset(mirna_ann, Species.Scientific.Name == "Homo sapiens")
mirna_name <- mirna_ann$Transcript.ID.Array.Design.
names(mirna_name) <- rownames(mirna_ann)

# Mature-miRNA-only filter (see "Analysis scope" above)
if (mature_only) {
  mature_ids <- rownames(mirna_ann)[mirna_ann$Sequence.Type == "miRNA"]
  mirna <- mirna[rownames(mirna) %in% mature_ids, ]
  cat("After mature-miRNA-only filter:", nrow(mirna), "probes\n")
}

# Rename rows: probeID_miRNAname
rownames(mirna) <- paste(rownames(mirna), mirna_name[rownames(mirna)], sep = "_")

# Rename columns: extract sample ID (4th underscore-delimited field)
colnames(mirna) <- sapply(strsplit(colnames(mirna), "_"), function(i) i[[4]])

# --- Intersect samples ---
comsubj <- intersect(rownames(clinical), colnames(mirna))
clinical <- clinical[comsubj, ]
clinical$Group <- factor(clinical$Group,
                          levels = c("Control", "Ambulatory", "Non-ambulatory"))
mirna <- mirna[, comsubj]

cat("Samples in common:", length(comsubj), "\n")
cat("Probes before filtering:", nrow(mirna), "\n")
cat("Sample alignment check:", all(rownames(clinical) == colnames(mirna)), "\n")

# --- Probe filtering (exact Rmd logic) ---
# 1. Keep probes with average expression > 3
mirna <- mirna[rowMeans(mirna) > 3, ]
cat("After rowMeans > 3 filter:", nrow(mirna), "probes\n")

# 2. Remove probes with NA annotation
mirna <- mirna[!(rownames(mirna) %in% grep("_NA", rownames(mirna), value = TRUE)), ]
cat("After removing _NA probes:", nrow(mirna), "probes\n")

# 3. Remove spike-in controls
mirna <- mirna[!(rownames(mirna) %in% grep("_spike", rownames(mirna), value = TRUE)), ]
cat("After removing _spike probes:", nrow(mirna), "probes\n")

# --- Save outputs ---
cat("Saving filtered data...\n")
write.csv(clinical, here("data", "clinical", "clinical.csv"))
write.csv(mirna, here("data", "mirna", "mirna.csv"))
saveRDS(mirna, here("results", "01_mirna_filtered.rds"))
saveRDS(clinical, here("results", "01_clinical.rds"))

cat("Done. Final dimensions:", nrow(mirna), "probes x", ncol(mirna), "samples\n")
