from rocker/tidyverse

ENV PATH="/home/miniconda3/bin:${PATH}"
ARG PATH="/home/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget libglpk-dev && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /home/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /home/miniconda3/ \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN Rscript -e "install.packages(c('xgboost', 'MCMCpack', 'reticulate', 'slider', 'vctrs', 'igraph', 'abind', 'ggh4x', 'fido'))"
RUN Rscript -e "install.packages('BiocManager'); BiocManager::install(c('phyloseq', 'DESeq2')); devtools::install_github('ruochenj/mbImpute/mbImpute R package');"

ENV HOME="/home/"
RUN chmod -R 777 /home/ /usr/local/