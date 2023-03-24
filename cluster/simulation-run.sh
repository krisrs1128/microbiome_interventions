
# docker run --user $(id -u):$(id -g) --rm=true -it -v $(pwd):/scratch -w /scratch 19e1c2fc5a7d /bin/bash


# install packages
tar -zxvf microbiome_interventions.tar.gz
cd microbiome_interventions
Rscript -e "purrr::map(c('mbtransfer', 'mdsine', 'tfPaper'), devtools::install)"
Rscript -e "mdsine::install_mdsine()"

# copy over data  
cp /staging/ksankaran/microbiome_interventions/tf_sim.tar.gz .
tar -zxvf tf_sim.tar.gz

# run the model configuration
PROCESS=$((PROCESS - 1))
for i in $(seq $((10 * PROCESS + 1)) $((10 * (PROCESS + 1)))); do
  export RUN=$(printf %03d $i)
  Rscript -e "rmarkdown::render('scripts/simulation_metrics.Rmd', params = list(data = 'tf_sim/sim_input_${RUN}.rda', run_id=${i}))"
done

mkdir result-${PROCESS}
mv result*rda result-${PROCESS}
tar -zcvf tf_sim_result-${PROCESS}.tar.gz result-${PROCESS}
cp tf_sim_result*tar.gz /staging/ksankaran/microbiome_interventions/