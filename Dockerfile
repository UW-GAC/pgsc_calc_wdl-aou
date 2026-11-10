# Software installation, no database files
FROM --platform=linux/amd64 condaforge/miniforge3:24.7.1-2

# build and run as root users since micromamba image has 'mambauser' set as the $USER
USER root

# set workdir to default for building; set to /data at the end
WORKDIR /

# install base dependencies; cleanup apt garbage
RUN apt-get update && apt-get install -y --no-install-recommends \
    bzip2 \
    ca-certificates \
    curl \
    git \
    gnupg2 \
    squashfs-tools \
    unzip \
    wget && \
    apt-get autoclean && \
    rm -rf /var/lib/apt /var/lib/dpkg /var/lib/cache /var/lib/log

# https://github.com/PGScatalog/pgsc_calc/discussions/302
RUN mamba create -y --name pgsc -c conda-forge -c bioconda -c defaults \
    python=3.10 \
    bioconda::plink2==2.00a3.3 \
    zstd \
    r-jsonlite \
    r-dplyr \
    r-tidyr \
    r-purrr \
    r-ggplot2 \
    r-DT \
    r-tibble \
    r-forcats \
    r-readr \
    quarto=1.4.550 \
    pip \
    nextflow=24.04.4 && \
    mamba clean -a -y

# set locale settings to UTF-8
# set the environment, put new hmas conda env in PATH by default
ENV PATH /opt/conda/envs/pgsc/bin:$PATH \
    LC_ALL=C.UTF-8

# set final working directory to /data
WORKDIR /data

# run test profile
# RUN nextflow run pgscatalog/pgsc_calc -r v2.0.0-alpha.5 -profile test,conda
