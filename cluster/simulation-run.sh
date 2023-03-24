
# docker run --user $(id -u):$(id -g) --rm=true -it -v $(pwd):/scratch -w /scratch 19e1c2fc5a7d /bin/bash


# install packages
Rscript -e "purrr::map(c('mbtransfer', 'mdsine', 'tfPaper'), devtools::install); mdsine::install_mdsine()"
Rscript -e "mdsine::install_mdsine()"
cd scripts

# run the model configuration
PROCESS=$((PROCESS - 1))
for i in $(seq $((10 * PROCESS + 1)) $((10 * (PROCESS + 1)))); do
  export RUN=$(printf %03d $i)
  cp /staging/groups/sankaran_group/microbiome_interventions/sim_input_${RUN}.rda .
  Rscript -e "rmarkdown::render('simulation_metrics.Rmd', params = list(data = 'sim_input_${RUN}.rda', run_id=${i}))"
done

