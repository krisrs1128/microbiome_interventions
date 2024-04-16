#!/bin/bash
# docker run --user $(id -u):$(id -g) --rm=true -it -v $(pwd):/scratch -w /scratch 19e1c2fc5a7d /bin/bash
  
tar -zxvf outputs.tar.gz

# copy and check if results already exist
test_output=outputs/forecasting-${process}.tar.gz
if [ -f $test_output ]; then
  cp $test_output .
  tar -zxvf forecasting-${process}.tar.gz || true
  file_num=$(ls -l forecasting-${process}/*.rda | wc -l) || true
else
  file_num=0
fi

if [[ $file_num -lt 1 ]]; then

  # install packages
  export batch_size=1
  tar -zxvf microbiome_interventions.tar.gz
  cd microbiome_interventions
  Rscript -e "install.packages('fido', repos='https://cloud.r-project.org')"
  Rscript -e "devtools::install_github('krisrs1128/microbiome_interventions/tfPaper')"
  Rscript -e "devtools::install_github('krisrs1128/microbiome_interventions/fido2')"
  Rscript -e "devtools::install_github('gathanei/xyz')"
  Rscript -e "devtools::install_github('krisrs1128/mbtransfer')"
  Rscript -e "devtools::install_github('krisrs1128/microbiome_interventions/mdsine')"
  Rscript -e "mdsine::install_mdsine()"

  # copy over data
  cp /staging/ksankaran/microbiome_interventions/tf_sim.tar.gz .
  tar -zxvf tf_sim.tar.gz

  # run the model configuration
  for i in $(seq $((batch_size * process + 1)) $((batch_size * (process + 1)))); do
    export RUN=$(printf %03d $i)
    Rscript -e "rmarkdown::render('scripts/forecasting_metrics.Rmd', params = list(data_dir = '../tf_sim/', run_id=${RUN}))"
  done

  mkdir forecasting-${process}
  mv scripts/*rda forecasting-${process}
  mv scripts/*html forecasting-${process}
  tar -zcvf forecasting-${process}.tar.gz forecasting-${process}
fi