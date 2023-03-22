
# docker run --user $(id -u):$(id -g) --rm=true -it -v $(pwd):/scratch -w /scratch 5713f1a635aa /bin/bash

Rscript -e "purrr::map(c('mbtransfer', 'mdsine', 'tfPaper'), devtools::install); mdsine::install_mdsine()"
Rscript -e "mdsine::install_mdsine()"

cd scripts
Rscript -e "rmarkdown::render('simulation_data.Rmd')"
Rscript -e "rmarkdown::render('simulation_metrics.Rmd')"