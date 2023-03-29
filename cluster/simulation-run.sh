
# docker run --user $(id -u):$(id -g) --rm=true -it -v $(pwd):/scratch -w /scratch 19e1c2fc5a7d /bin/bash

# install packages
export batch_size=10
tar -zxvf microbiome_interventions.tar.gz
cd microbiome_interventions
Rscript -e "purrr::map(c('mbtransfer', 'mdsine', 'tfPaper'), devtools::install)"
Rscript -e "mdsine::install_mdsine()"

# copy over data  
cp /staging/ksankaran/microbiome_interventions/tf_sim.tar.gz .
tar -zxvf tf_sim.tar.gz

# run the model configuration
for i in $(seq $((batch_size * process + 1)) $((batch_size * (process + 1)))); do
  export RUN=$(printf %03d $i)
  Rscript -e "rmarkdown::render('scripts/simulation_metrics.Rmd', params = list(data = 'tf_sim/sim_input_${RUN}.rda', run_id=${i}))"
done

mkdir result-${process}
mv result*rda result-${process}
tar -zcvf tf_sim_result-${process}.tar.gz result-${process}
cp tf_sim_result*tar.gz /staging/ksankaran/microbiome_interventions/
