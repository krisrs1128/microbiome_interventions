from rocker/tidyverse

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

RUN Rscript -e "install.packages('BiocManager'); BiocManager::install(c('phyloseq', 'DESeq2')); devtools::install_github('ruochenj/mbImpute/mbImpute R package')"
RUN chmod -R 777 /usr/local/lib
