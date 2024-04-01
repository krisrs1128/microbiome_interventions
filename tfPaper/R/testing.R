#' @export
testing_metrics <- function(R, truth) {
  S <- intersect(truth, R)
  V <- setdiff(R, truth)

  list(
    fdp = length(V) / max(length(R), 1),
    power = length(S) / length(truth)
  )
}

#' @importFrom dplyr bind_rows
#' @export
lagged_testing_metrics <- function(R, nonnull_taxa) {
  testing_results <- list()
  for (l in seq_along(R)) {
    testing_results[[l]] <- list()
    testing_results[[l]]$marginal <- testing_metrics(R[[l]], nonnull_taxa[[l]]$marginal)
    testing_results[[l]]$interaction <- testing_metrics(R[[l]], nonnull_taxa[[l]]$interaction)
    testing_results[[l]]$interaction$fdp <- NA

    testing_results[[l]] <- testing_results[[l]] |>
      bind_rows(.id = "type")
  }

  bind_rows(testing_results, .id = "lag")
}
