#!/bin/bash
# docker run --user $(id -u):$(id -g) --rm=true -it -v $(pwd):/scratch -w /scratch 19e1c2fc5a7d /bin/bash

# install packages
export batch_size=1
tar -zxvf microbiome_interventions.tar.gz
cd microbiome_interventions
Rscript -e "purrr::map(c('mbtransfer', 'mdsine', 'tfPaper', 'fido2'), devtools::install)"
Rscript -e "mdsine::install_mdsine()"

# copy over data  
cp /staging/ksankaran/microbiome_interventions/tf_sim.tar.gz .
tar -zxvf tf_sim.tar.gz

# used to rerun high memory jobs
export run_ids=('271' '275' '279' '283' '287' '291' '295' '299' '303' '307' '311' '315' '319' '323' '327' '331' '335' '339' '343' '347' '351' '355' '359' '363' '367' '371' '375')
process=${run_ids[process]}

# run the model configuration
for i in $(seq $((batch_size * process + 1)) $((batch_size * (process + 1)))); do
  export RUN=$(printf %03d $i)
  Rscript -e "rmarkdown::render('scripts/forecasting_metrics.Rmd', params = list(data_dir = '../tf_sim/', run_id=${RUN}))"
done

mkdir forecasting-${process}
mv scripts/*rda forecasting-${process}
mv scripts/*html forecasting-${process}
tar -zcvf forecasting-${process}.tar.gz forecasting-${process}
cp forecasting-*tar.gz /staging/ksankaran/microbiome_interventions/
