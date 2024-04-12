
# docker run --user $(id -u):$(id -g) --rm=true -it -v $(pwd):/scratch -w /scratch 19e1c2fc5a7d /bin/bash

# copy and check if results already exist
cp /staging/ksankaran/microbiome_interventions/inference-${process}.* . || true
tar -zxvf inference-${process}.tar.gz || true

file_num=$(ls -l inference-${process}/*.rda | wc -l) || true
if [[ $file_num -lt 1 ]]; then

  # install packages
  export batch_size=1
  tar -zxvf microbiome_interventions.tar.gz
  cd microbiome_interventions
  Rscript -e "devtools::install('tfPaper')"
  Rscript -e "devtools::install_github('gathanei/xyz')"
  Rscript -e "devtools::install_github('krisrs1128/mbtransfer')"

  # copy over data  
  cp /staging/ksankaran/microbiome_interventions/tf_sim.tar.gz .
  tar -zxvf tf_sim.tar.gz

  #process=$((process + 54))

  # run the model configuration
  for i in $(seq $((batch_size * process + 1)) $((batch_size * (process + 1)))); do
    export RUN=$(printf %03d $i)
    Rscript -e "rmarkdown::render('scripts/inference_metrics.Rmd', params = list(data_dir = '../tf_sim/', run_id=${RUN}))"
  done

  mkdir inference-${process}
  mv scripts/*rda inference-${process}
  mv scripts/*html inference-${process}
  tar -zcvf inference_result-${process}.tar.gz inference-${process}
  cp inference_result*tar.gz /staging/ksankaran/microbiome_interventions/

fi;