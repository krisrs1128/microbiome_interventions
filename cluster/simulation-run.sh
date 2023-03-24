
# docker run --user $(id -u):$(id -g) --rm=true -it -v $(pwd):/scratch -w /scratch 19e1c2fc5a7d /bin/bash

# copy the data over
export RUN=$(printf %03d $PROCESS)
cp /staging/groups/sankaran_group/microbiome_interventions/sim_input_${RUN}.rda .

# install packages
Rscript -e "purrr::map(c('mbtransfer', 'mdsine', 'tfPaper'), devtools::install); mdsine::install_mdsine()"
Rscript -e "mdsine::install_mdsine()"

# run the model configuration
Rscript -e "rmarkdown::render('scripts/simulation_metrics.Rmd', params = list(data = 'sim_input_${RUN}.rda', run_id=${PROCESS}))"