FROM rocker/rstudio:4.4.2

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libsqlite3-dev \
    zlib1g-dev \
    libglpk-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Set CRAN mirror and Bioconductor version
ENV CRAN=https://cran.r-project.org
ENV BIOCONDUCTOR_VERSION=3.20

# Install BiocManager and set version
RUN R -e "install.packages('BiocManager', repos='https://cran.r-project.org'); \
  BiocManager::install(version='3.20', ask=FALSE, update=FALSE)"

# Bioconductor packages
# Install preprocessCore first with --disable-threading because oligo depends on it.
RUN R -e "BiocManager::install('preprocessCore', \
  configure.args='--disable-threading', ask=FALSE, update=FALSE, force=TRUE); \
  if (!requireNamespace('preprocessCore', quietly=TRUE)) stop('preprocessCore install failed')"

RUN R -e "BiocManager::install(c('oligo', 'pd.mirna.4.0'), \
  ask=FALSE, update=FALSE); \
  missing <- c('oligo', 'pd.mirna.4.0')[!vapply(c('oligo', 'pd.mirna.4.0'), requireNamespace, logical(1), quietly=TRUE)]; \
  if (length(missing)) stop('Missing Bioconductor package(s): ', paste(missing, collapse=', '))"

RUN R -e "BiocManager::install(c('limma', 'ComplexHeatmap'), \
  ask=FALSE, update=FALSE)"

# CRAN packages
RUN R -e "install.packages(c( \
  'tidyverse', 'ggpubr', 'cowplot', 'ggrepel', 'RColorBrewer', \
  'enrichR', 'RSQLite', 'here', 'circlize', 'igraph', \
  'tidygraph', 'ggraph' \
  ))"

# Create working directory
WORKDIR /home/rstudio

# Copy .Rprofile
COPY .Rprofile /home/rstudio/.Rprofile

# Set default project
ENV RSTUDIO_DEFAULT_WORKING_DIRECTORY=/home/rstudio
