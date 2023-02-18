
#' importFrom cmdstanr cmdstan_model
#' @export
train <- function(ts_inter, method = "zeros", hyper = list()) {
  if (method == "zeros") {
    result <- new("transfer_model", parameters = list(), method = method)
  }
  if (method == "gaussian_latent") {
    model <- cmdstan_model(system.file("gaussian_latent.stan", package = "microTF"))
    data_list <- list(
      x = ts_inter[[1]]@values,
      interventions = ts_inter[[1]]@interventions,
      P = hyper$P,
      Q = hyper$Q,
      K = hyper$K,
      D = nrow(ts_inter[[1]]@interventions),
      N = ncol(ts_inter[[1]]@values),
      V = nrow(ts_inter[[1]]@values)
    )
    fit <- model$variational(data_list)
    result <- new("transfer_model", parameters = fit, method = method)
  }
  
  result
}

zeros_predict <- function(newdata) {
  for (i in seq_along(newdata)) {
    preds <- matrix(0, nrow(newdata[[i]]), ncol(newdata[[i]]@interventions) - ncol(newdata[[i]]))
    values(newdata[[i]]) <- cbind(values(newdata[[i]]), preds)
  }
  newdata
}
