import mdsine2 as md2
import pathlib
from mdsine2.names import STRNAMES

def set_params(model="NegBin", **kwargs):
  params = {"seed": 0, "burnin": 0, "n_samples": 20, "checkpoint": 10, "basepath": "."}
  params.update(**kwargs)

  if model == "NegBin":
    result = md2.config.NegBinConfig(**params)
  else:
    result = md2.config.MDSINE2ModelConfig(**params)
  return result


def md_data(paths):
  return md2.dataset.parse(
    name="mdsine-model", 
    taxonomy=paths["taxonomy"],
    reads=paths["reads"],
    qpcr=paths["qpcr"],
    metadata=paths["metadata"],
    perturbations=paths["perturbations"]
  )
  
def mdsine(dataset, **kwargs):
  params = set_params(model="NegBin", **kwargs)
  basepath = pathlib.Path(params.__dict__["OUTPUT_BASEPATH"])
  basepath.mkdir(exist_ok=True, parents=True)
  
  # setup the negative binomial initial run
  mcmc_negbin = md2.negbin.build_graph(
      params=params,
      graph_name=dataset.name,
      subjset=dataset
  )
  mcmc_negbin.run()

  # setup parameters for the full run
  a0 = md2.summary(mcmc_negbin.graph[STRNAMES.NEGBIN_A0])["mean"]
  a1 = md2.summary(mcmc_negbin.graph[STRNAMES.NEGBIN_A1])["mean"]
  params = set_params(model="MCMC", negbin_a0=a0, negbin_a1=a1, **kwargs)
  mcmc = md2.initialize_graph(params=params, graph_name=dataset.name, subjset=dataset)
  return md2.run_graph(mcmc)
