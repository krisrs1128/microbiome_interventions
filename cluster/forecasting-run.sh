
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
export run_ids=('271' '272' '275' '276' '279' '280' '283' '284' '287' '288' '291' '292' '295' '296' '299' '300' '303' '304' '307' '308' '311' '312' '315' '316' '319' '320' '323' '324' '327' '328' '331' '332' '335' '336' '339' '340' '343' '344' '347' '348' '351' '352' '355' '356' '359' '360' '363' '364' '367' '368' '371' '372' '375' '376')
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
