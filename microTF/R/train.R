
#' @importFrom xgboost xgboost
#' @importFrom glmnet glmnet
#' @export
train <- function(ts_inter, method = "zeros", hyper = list()) {
  if (method == "zeros") {
    result <- new("transfer_model", parameters = list(), method = method)
  } else if (method == "gaussian_latent") {

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
    
  } else if (method == "mdsine") {
    
    fit <- mdsine(ts, hyper$taxonomy)
    result <- new("transfer_model", parameters = fit, method = method)
    
  } else if(method == "gbm") {
    train_data <- patchify_df(ts_inter)
    fit <- list()
    for (j in seq_along(train_data$y)) {
      fit[[j]] <- xgboost(data = train_data$x, label = train_data$y[[j]], nrounds=50, verbose=0, nthread=4)
    }
    result <- new("transfer_model", parameters = fit, method = method)

  } else if (method == "lasso") {
    train_data <- patchify_df(ts_inter)
    fit <- list()
    for (j in seq_along(train_data$y)) {
      fit[[j]] <- glmnet(as.matrix(train_data$x), train_data$y[[j]])
    }
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
