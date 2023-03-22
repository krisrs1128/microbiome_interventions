
# docker run --user $(id -u):$(id -g) --rm=true -it -v $(pwd):/scratch -w /scratch (my image ID) /bin/bash

tar -zxvf microbiome_interventions.tar.gz
cd microbiome_interventions
Rscript -e "install.packages('BiocManager'); BiocManager::install(c('phyloseq', 'DESeq2')); devtools::install_github('ruochenj/mbImpute/mbImpute R package')"
Rscript -e "purrr::map(c('tfPaper', 'mbtransfer', 'mdsine'), devtools::install); mdsine::install_mdsine()"
Rscript -e "mdsine::install_mdsine()"

cd scripts
Rscript -e "rmarkdown::render('simulation_data.Rmd')"
Rscript -e "rmarkdown::render('simulation_metrics.Rmd')"