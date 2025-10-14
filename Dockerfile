# Pin to the latest R 4.4 patch for stability
FROM rocker/rstudio:4.4.1

LABEL maintainer="kdreval@sfu.ca"
LABEL description="Docker container for GAMBLR.predict (pinned & low-RAM build)"

# ---- Versions you control (tags or SHAs) ----
ARG GAMBLR_DATA_REF=v1.3.1
ARG GAMBLR_HELPERS_REF=v1.3.1
ARG GAMBLR_PREDICT_REF=v1.2

# ---- Keep RAM use low during installs ----
ENV PAK_NUM_WORKERS=1 \
    MAKEFLAGS=-j1 \
    R_INSTALL_STAGED=false

# ---- System deps (single layer, no recommends, then clean) ----
USER root
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    libxml2-dev libssl-dev libcurl4-openssl-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    unixodbc-dev cmake \
 && rm -rf /var/lib/apt/lists/*

# ---- Use pak for deterministic installs ----
ARG CRAN_MIRROR=https://cran.r-project.org
ENV RENV_CONFIG_REPOS_OVERRIDE=${CRAN_MIRROR}
RUN R -q -e "install.packages('pak', repos='${CRAN_MIRROR}'); pak::pak_setup()"
RUN R -q -e "install.packages('BiocManager', repos='${CRAN_MIRROR}')"

# ---- CRAN packages ----
RUN R -q -e "pak::pkg_install(c( \
  'caret','circlize','dplyr','DT','FNN','ggalluvial','ggrepel','ggside', \
  'grid','mclust','purrr','randomForest','readr','shinybusy','shinyjs', \
  'tibble','tidyr','tidyselect','plotly','ggpubr','readxl','uwot' \
), ask = FALSE)"

# ---- Bioconductor packages (from Bioc, not GitHub) ----
RUN R -q -e "BiocManager::install(c('IRanges','ComplexHeatmap'), ask=FALSE, update=FALSE)"

# Install remotes (small, fast), then install GAMBLR.data with low-memory flags
RUN R -q -e "install.packages('remotes', repos='https://cran.r-project.org')"
RUN R -q -e " \
  ref <- sprintf('morinlab/GAMBLR.data@%s', Sys.getenv('GAMBLR_DATA_REF', 'v1.3.1')); \
  remotes::install_github(ref, upgrade = 'never', \
    INSTALL_opts = c('--no-multiarch','--no-byte-compile','--no-build-vignettes','--no-test-load'), \
    Ncpus = 1); \
  stopifnot(requireNamespace('GAMBLR.data', quietly = TRUE)); \
  cat('GAMBLR.data ', as.character(packageVersion('GAMBLR.data')), '\n') \
"

# Pin uwot if you need a specific GitHub version (CRAN uwot is usually fine)
# RUN R -q -e "pak::pkg_install('jlmelville/uwot@v0.2.3', ask=FALSE)"

RUN R -q -e "pak::pkg_install(sprintf('morinlab/GAMBLR.helpers@%s', Sys.getenv('GAMBLR_HELPERS_REF', '${GAMBLR_HELPERS_REF}')), ask=FALSE)"
RUN R -q -e "pak::pkg_install(sprintf('morinlab/GAMBLR.predict@%s', Sys.getenv('GAMBLR_PREDICT_REF', '${GAMBLR_PREDICT_REF}')), ask=FALSE)"

# ---- Smoke tests: fail early if anything is missing ----
RUN R -q -e " \
  stopifnot(requireNamespace('GAMBLR.data',     quietly=TRUE)); \
  stopifnot(requireNamespace('GAMBLR.helpers',  quietly=TRUE)); \
  stopifnot(requireNamespace('GAMBLR.predict',  quietly=TRUE)); \
  cat('Installed versions:\\n'); \
  for(p in c('GAMBLR.data','GAMBLR.helpers','GAMBLR.predict')) \
    cat('  ', p, as.character(packageVersion(p)), '\\n') \
"
RUN R -q -e "pak::pkg_install('jlmelville/uwot@v0.2.3', ask=FALSE)"
EXPOSE 8787 3838

HEALTHCHECK CMD R -q -e "cat('ok')" || exit 1


