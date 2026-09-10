# DYSF miRNA Biomarker — Reproducible Docker Analysis

Reproduce the entire DYSF miRNA analysis from CEL files through all 6 manuscript figures in a self-contained RStudio environment.

## Quick Start

```bash
# 1. Download CEL files from GSE346695 and place them in the cel_files/ directory
cp /path/to/CEL/files/*.CEL cel_files/

# 2a. Build and launch RStudio
docker compose up -d --build
# or 2b. Build and launch RStudio without cache
docker compose build --no-cache && docker compose up -d

# 3. Open RStudio in your browser
#    URL: http://localhost:8787
#    Username: rstudio
#    Password: dysf2025

# 4. Run scripts in order (in RStudio console or terminal)
Rscript src/00_cel_to_normalized.R    # CEL → normalized data (one-time)
Rscript src/00b_qc_cel.R              # CEL quality-control report
Rscript src/01_data_import.R          # Load + filter → 666 probes
Rscript src/02_FigS1_pca_heatmap.R    # Figure S1: PCA × clinical heatmap
Rscript src/03_Fig1_fdr_volcano.R     # Figure 1: FDR rank + volcano plots
Rscript src/04_Fig2A_mir4532boxplot.R # Figure 2A: miR-4532 boxplot
Rscript src/05_Fig2B_cbc_heatmap.R    # Figure 2B: CBC correlation heatmap
Rscript src/06_Fig2C_monocyte_scatter.R # Figure 2C: Monocyte scatter
Rscript src/07_Fig2D_enrichr_network.R  # Figure 2D: enrichR network/dotplot
```

## Script Descriptions

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `00_cel_to_normalized.R` | `cel_files/*.CEL` | `data/mirna/normalized_data_celderived.TXT` | RMA normalization via oligo (pd.mirna.4.0, checkType=FALSE) |
| `00b_qc_cel.R` | `cel_files/*.CEL` | `results/qc/cel_qc_report.pdf`, QC PNGs and CSVs | CEL-level quality-control diagnostics |
| `01_data_import.R` | `normalized_data_celderived.TXT`, clinical/CBC CSVs | `results/01_mirna_filtered.rds`, `results/01_clinical.rds` | Load data, annotate probes, filter to 666 probes |
| `02_FigS1_pca_heatmap.R` | RDS from step 1 | `results/FigS1_pca_heatmap.png` | PCA (rank=5) × clinical variable heatmap (ANOVA p-values) |
| `03_Fig1_fdr_volcano.R` | RDS from step 1 | `results/Fig1_fdr_volcano.png` | limma DE (3 contrasts), FDR rank plot + volcano |
| `04_Fig2A_mir4532boxplot.R` | RDS from steps 1–2 | `results/Fig2A_mir4532boxplot.png` | Two-stage candidate selection, miR-4532 boxplot |
| `05_Fig2B_cbc_heatmap.R` | RDS from steps 1–2 | `results/Fig2B_cbc_heatmap.pdf` | Per-group Spearman CBC × miR-4532 heatmap |
| `06_Fig2C_monocyte_scatter.R` | RDS from step 1 | `results/Fig2C_monocyte_scatter.pdf` | Monocyte × miR-4532 Pearson scatter per group |
| `07_Fig2D_enrichr_network.R` | `mir2gene.sqlite` | `results/Fig2D_enrichr_network.pdf` or `results/Fig2D_enrichr_dotplot.png` | Target enrichment (enrichR), network or dotplot fallback |

## Notes

- **Script 00 is optional.** The downstream scripts use the pre-existing `data/mirna/normalized_data_celderived.TXT`. Script 00 is provided to reproduce the CEL-to-normalization step from scratch.
- **CEL files are not bundled.** Mount them separately via the `cel_files/` volume.
- **enrichR requires internet access.** The container must be able to reach the enrichR API. If no terms survive FDR < 0.1, a dotplot of top nominal terms is generated as fallback.
- **R version:** R 4.4.2 / Bioconductor 3.20, pinned for reproducibility.
- **preprocessCore** is installed with `--disable-threading` to prevent `pthread_create() error 22` on some environments.

## Directory Structure

```
dysf_mirna/
├── .gitignore
├── .Rprofile
├── .here
├── Dockerfile
├── README.md
├── docker-compose.yml
├── cel_files/          # Raw CEL files (74 arrays; not tracked by Git)
├── data/
│   ├── annotations/
│   │   └── mir2gene.sqlite
│   ├── cbc/
│   │   └── md_cbc.csv
│   ├── clinical/
│   │   ├── Samples sent to Scripps_2015.04.07_age_sex.csv
│   │   └── clinical.csv
│   └── mirna/
│       ├── miRNA-4_0-st-v1.annotations.20160922.csv
│       ├── mirna.csv
│       └── normalized_data_celderived.TXT
├── geo/
│   └── GA_affy_DYSF.xlsx
├── manuscript/
│   ├── old/
│   ├── tgrewal_dysf_mirna_manuscript_2026.docx
│   └── tgrewal_dysf_mirna_supplement_2026.docx
├── results/            # Figures, QC reports, and RDS intermediates
│   └── qc/
└── src/
│   ├── utils.R
│   ├── 00_cel_to_normalized.R
│   ├── 00b_qc_cel.R
│   ├── 01_data_import.R
│   ├── 02_FigS1_pca_heatmap.R
│   ├── 03_Fig1_fdr_volcano.R
│   ├── 04_Fig2A_mir4532boxplot.R
│   ├── 05_Fig2B_cbc_heatmap.R
│   ├── 06_Fig2C_monocyte_scatter.R
│   └── 07_Fig2D_enrichr_network.R
```
