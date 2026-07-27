# Dockerfile for R-workflows

# Use rocker/tidyverse as base image for common bioinformatics R packages
FROM rocker/tidyverse:4.2.3

# Install system dependencies commonly needed for Bioconductor and other packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    && rm -rf /var/lib/apt/lists/*

# Install remotes and renv for reproducible environments
RUN R -e "install.packages(c('remotes','renv'), repos='https://cloud.r-project.org')"

# Create workdir
WORKDIR /work

# Copy project files (use a small context when building to avoid copying large data)
COPY . /work

# Default command prints usage. To run a script use: docker run --rm -v $(pwd):/work -w /work r-workflows:latest Rscript r-00_rawdata.R
CMD ["/bin/bash"]
