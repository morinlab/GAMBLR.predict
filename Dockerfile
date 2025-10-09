# Use a recent R version from rocker
FROM rocker/rstudio:4.4.0

# Add labels
LABEL maintainer="kdreval@sfu.ca"
LABEL description="Docker container for GAMBLR.predict"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    unixodbc-dev

# Install remotes package first to install other dependencies
RUN R -e "install.packages('remotes', dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install CRAN packages
RUN R -e "remotes::install_cran(c('caret', 'circlize', 'dplyr', 'DT', 'FNN', 'ggalluvial', 'ggrepel', 'ggside', 'grid', 'mclust', 'purrr', 'randomForest', 'readr', 'shinybusy', 'shinyjs', 'tibble', 'tidyr', 'tidyselect', 'uwot'))"

# Install Bioconductor packages
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('IRanges',ask=FALSE)"

# Install GitHub packages

RUN R -e "remotes::install_github('jokergoo/ComplexHeatmap')"
RUN R -e "remotes::install_github('morinlab/GAMBLR.data')"
RUN R -e "remotes::install_github('morinlab/GAMBLR.helpers')"
RUN apt-get install -y cmake 
RUN R -e "install.packages('ggpubr',dependencies=TRUE)"
RUN R -e "install.packages('plotly',dependencies=TRUE)"
RUN R -e "remotes::install_github('jlmelville/uwot')"
RUN R -e "remotes::install_github('morinlab/GAMBLR.predict')"
