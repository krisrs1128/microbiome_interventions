#!/bin/bash
# docker run --user $(id -u):$(id -g) --rm=true -it -v $(pwd):/scratch -w /scratch 19e1c2fc5a7d /bin/bash
  
# install packages
tar -zxvf microbiome_interventions.tar.gz
tar -zxvf tf_sim.tar.gz
cd microbiome_interventions
Rscript -e "install.packages('fido', repos='https://cloud.r-project.org')"
Rscript -e "devtools::install_github('krisrs1128/microbiome_interventions/tfPaper')"
Rscript -e "devtools::install_github('krisrs1128/microbiome_interventions/fido2')"
Rscript -e "devtools::install_github('gathanei/xyz')"
Rscript -e "devtools::install_github('krisrs1128/mbtransfer')"
Rscript -e "devtools::install_github('krisrs1128/microbiome_interventions/mdsine')"
Rscript -e "mdsine::install_mdsine()"

# run the model configuration
Rscript -e "source(knitr::purl('scripts/compare_model_trajectories.Rmd'))"

mkdir comparison-run
mv *rda comparison-run
mv *html comparison-run
mv *rds comparison-run
tar -zcvf comparison-run.tar.gz comparison-run
cp comparison-run.tar.gz ../
