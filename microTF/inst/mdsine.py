import mdsine2 as md2

def default_params():
  params = {}
  params["nb"] = md2.config.NegBinConfig(seed=0, burnin=100, n_samples=200, checkpoint=50, basepath="output")
  params["full"] = md2.config.MDSINE2ModelConfig(seed=0, burnin=50, n_samples=100, checkpoint=25, basepath="output")
  return params


def fit_mdsine(params, dataset, basepath="output"):
  basepath = pathlib.Path(basepath)
  basepath.mkdir(exist_ok=True, parents=True)
  
  # setup the negative binomial initial run
  mcmc_negbin = md2.negbin.build_graph(
      params=params["nb"],
      graph_name=dataset.name,
      subjset=dataset
  )
  mcmc_negbin.run()
  
  # setup parameters for the full run
  a0 = md2.summary(mcmc_negbin.graph[STRNAMES.NEGBIN_A0])["mean"]
  a1 = md2.summary(mcmc_negbin.graph[STRNAMES.NEGBIN_A1])["mean"]
  params["full"].update({
    basepath=str(basepath),
    negbin_a0=a0, 
    negbin_a1=a1
  })

  mcmc = md2.initialize_graph(params=params, graph_name=dataset.name, subjset=dataset)
  return md2.run_graph(mcmc)
