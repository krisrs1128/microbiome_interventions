
#' @importFrom dplyr mutate filter
#' @importFrom glue glue
#' @importFrom tidyr separate
#' @export
method_configurations <- function(data_paths) {
  expand.grid(
    method_hyper = c("mbtransfer-1", "mbtransfer-2", "mdsine"),
    normalization = c("none", "DESeq2", "mbImpute"),
    data_path = data_paths
  ) |>
    separate(method_hyper, c("method", "hyper"), convert = TRUE) |>
    mutate(
      hyper = list(list(P = 2, Q = 2), list(P = 4, Q = 4))[hyper],
      output_path = glue("result-{str_pad(row_number(), 2, 'left', '0')}.rda")
    ) |>
    filter(!(normalization == "mbImpute" & method == "mdsine"))
}

#' @importFrom dplyr mutate
#' @importFrom glue glue
#' @export
varying_parameters <- function(output_root="~/Downloads/") {
  expand.grid(
    n_subject = c(20, 40),
    n_time = c(20, 40),
    prop_nonnull = c(0.05, 0.2),
    signal_B = c(0.5, 1, 2),
    replicate = seq_len(3)
  ) |>
    mutate(output_path = glue("{output_root}/sim_input_{str_pad(row_number(), 3, 'left', '0')}.rda"))
}

#' @export
fixed_parameters <- function() {
  list(
    n_perturb = 1,
    n_covariates = 3, 
    n_latent = 2,
    normalizer_A = 0.8,
    normalizer_C = 0.8,
    sparsity_A = 0.8,
    n_lag = 4,
    n_taxa = 250
  )
}

