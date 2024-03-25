
#' @importFrom dplyr mutate filter
#' @importFrom glue glue
#' @importFrom tidyr separate
#' @export
method_configurations <- function(data_paths) {
  expand.grid(
    method_hyper = c("mbtransfer-1", "mbtransfer-2", "mdsine", "fido-3", "fido-4"),
    normalization = c("none", "DESeq2", "DESeq2-asinh"),
    data_path = data_paths
  ) |>
    separate(method_hyper, c("method", "hyper"), convert = TRUE) |>
    mutate(
      hyper = list(list(P = 2, Q = 2), list(P = 4, Q = 4), list(sigma = 1, rho = 1), list(sigma = .5, rho = .5))[hyper],
      output_path = glue("result-{str_pad(row_number(), 3, 'left', '0')}.rda")
    )
}

#' @importFrom dplyr mutate
#' @importFrom glue glue
#' @export
data_parameters <- function(output_root="~/Downloads/") {
  expand.grid(
    n_subject = 50,
    n_time = 30,
    prop_nonnull = c(0.1, 0.2, 0.4),
    signal_B = c(0.25, 0.5, 1),
    n_taxa = c(100, 200, 400),
    phylo_alpha = c(0.1, 1, 10),
    baseline_lambda = c(0.1, 1, 10)
  ) |>
    mutate(output_path = glue("{output_root}/sim_input_{str_pad(row_number(), 3, 'left', '0')}.rda"))
}

#' @export
fixed_parameters <- function() {
  list(
    n_perturb = 1,
    n_covariates = 1, 
    prop_interaction = 0.4,
    n_latent = 4,
    normalizer_A = 0.5,
    normalizer_C = 0.2,
    sparsity_A = 0.9,
    signal_C = 0.5,
    n_lag = 3
  )
}
