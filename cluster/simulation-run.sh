
# docker run --user $(id -u):$(id -g) --rm=true -it -v $(pwd):/scratch -w /scratch 19e1c2fc5a7d /bin/bash


# install packages
Rscript -e "purrr::map(c('mbtransfer', 'mdsine', 'tfPaper'), devtools::install); mdsine::install_mdsine()"
Rscript -e "mdsine::install_mdsine()"
cd scripts

# copy over data  
cp /staging/ksankaran/microbiome_interventions/tf_sim.tar.gz .
tar -zxvf tf_sim.tar.gz

# run the model configuration
PROCESS=$((PROCESS - 1))
for i in $(seq $((10 * PROCESS + 1)) $((10 * (PROCESS + 1)))); do
  export RUN=$(printf %03d $i)
  Rscript -e "rmarkdown::render('simulation_metrics.Rmd', params = list(data = 'sim_input_${RUN}.rda', run_id=${i}))"
done

mkdir result-${PROCESS}
mv result* result-${PROCESS}
tar -zcvf tf_sim_result-${PROCESS}.tar.gz result-${PROCESS}
cp tf_sim_result*tar.gz /staging/ksankaran/microbiome_interventions/