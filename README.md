# DYSF miRNA Biomarker — Reproducible Docker Analysis

Reproduce the entire DYSF miRNA analysis from CEL files through all 6 manuscript figures in a self-contained RStudio environment.

## Quick Start

```bash
# 1. Place your CEL files in the cel_files/ directory
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
Rscript src/01_data_import.R          # Load + filter → 666 probes
Rscript src/02_figure1_pca_heatmap.R  # Figure 1: PCA × clinical heatmap
Rscript src/03_figure2_volcano.R      # Figure 2: FDR rank + volcano plots
Rscript src/04_figure3_mir4532_boxplot.R  # Figure 3: miR-4532 boxplot
Rscript src/05_figure4_cbc_heatmap.R  # Figure 4: CBC correlation heatmap
Rscript src/06_figure5_monocyte_scatter.R # Figure 5: Monocyte scatter
Rscript src/07_figure6_enrichr_network.R  # Figure 6: enrichR network/dotplot
```

## Script Descriptions

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `00_cel_to_normalized.R` | `cel_files/*.CEL` | `data/mirna/normalized_data.TXT` | RMA normalization via oligo (pd.mirna.4.0, checkType=FALSE) |
| `01_data_import.R` | `normalized_data.TXT`, clinical/CBC CSVs | `results/01_mirna_filtered.rds`, `results/01_clinical.rds` | Load data, annotate probes, filter to 666 probes |
| `02_figure1_pca_heatmap.R` | RDS from step 1 | `results/fig1_pca_heatmap.png` | PCA (rank=5) × clinical variable heatmap (ANOVA p-values) |
| `03_figure2_volcano.R` | RDS from step 1 | `results/fig2_volcano.png` | limma DE (3 contrasts), FDR rank plot + volcano |
| `04_figure3_mir4532_boxplot.R` | RDS from steps 1–2 | `results/fig3_mir4532_boxplot.png` | Two-stage candidate selection, miR-4532 boxplot |
| `05_figure4_cbc_heatmap.R` | RDS from steps 1–2 | `results/fig4_cbc_heatmap.pdf` | Per-group Spearman CBC × miR-4532 heatmap |
| `06_figure5_monocyte_scatter.R` | RDS from step 1 | `results/fig5_monocyte_scatter.pdf` | Monocyte × miR-4532 Pearson scatter per group |
| `07_figure6_enrichr_network.R` | `mir2gene.sqlite` | `results/fig6_enrichr_network.pdf` or `fig6_enrichr_dotplot.png` | Target enrichment (enrichR), network or dotplot fallback |

## Notes

- **Script 00 is optional.** The downstream scripts use the pre-existing `data/mirna/normalized_data.TXT`. Script 00 is provided to reproduce the CEL-to-normalization step from scratch.
- **CEL files are not bundled.** Mount them separately via the `cel_files/` volume.
- **enrichR requires internet access.** The container must be able to reach the enrichR API. If no terms survive FDR < 0.1, a dotplot of top nominal terms is generated as fallback.
- **R version:** R 4.4.2 / Bioconductor 3.20, pinned for reproducibility.
- **preprocessCore** is installed with `--disable-threading` to prevent `pthread_create() error 22` on some environments.

## Directory Structure

```
dysf_mirna_docker/
├── Dockerfile
├── docker-compose.yml
├── .Rprofile
├── .here
├── README.md
├── data/
│   ├── mirna/
│   │   ├── normalized_data.TXT
│   │   └── miRNA-4_0-st-v1.annotations.20160922.csv
│   ├── clinical/
│   │   └── Samples sent to Scripps_2015.04.07_age_sex.csv
│   ├── cbc/
│   │   └── md_cbc.csv
│   └── annotations/
│       └── mir2gene.sqlite
├── cel_files/          # YOU PROVIDE: 74 CEL files here
├── src/
│   ├── utils.R
│   ├── 00_cel_to_normalized.R
│   ├── 01_data_import.R
│   ├── 02_figure1_pca_heatmap.R
│   ├── 03_figure2_volcano.R
│   ├── 04_figure3_mir4532_boxplot.R
│   ├── 05_figure4_cbc_heatmap.R
│   ├── 06_figure5_monocyte_scatter.R
│   └── 07_figure6_enrichr_network.R
└── results/           # Output directory (figures + RDS intermediates)
```
