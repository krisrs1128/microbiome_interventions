#!/bin/bash

tar -zxvf microbiome_interventions.tar.gz
cd microbiome_interventions
Rscript -e "install.packages('fido', repos='https://cloud.r-project.org')"
Rscript -e "devtools::install_github('krisrs1128/microbiome_interventions/tfPaper')"
Rscript -e "devtools::install_github('krisrs1128/microbiome_interventions/fido2')"
Rscript -e "devtools::install_github('gathanei/xyz')"
Rscript -e "install.packages('prettyunits')"
Rscript -e "devtools::install_github('krisrs1128/mbtransfer')"
Rscript -e "devtools::install_github('krisrs1128/microbiome_interventions/mdsine')"
Rscript -e "mdsine::install_mdsine()"

Rscript -e "rmarkdown::render('case_studies/diet.Rmd')"
cp case_studies/*.rda /staging/ksankaran/