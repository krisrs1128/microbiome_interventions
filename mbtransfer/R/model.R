
#' @importFrom xgboost xgboost
#' @export
mbtransfer <- function(ts_inter, P = 1, Q = 1, verbose = 0, nrounds = 200, ...) {
  train_data <- patchify_df(ts_inter, P, Q)
  fit <- list()
  for (j in seq_along(train_data$y)) {
    fit[[j]] <- xgboost(data = train_data$x, label = train_data$y[[j]], nrounds = nrounds, verbose = verbose, ...)
  }
  
  hyper <- list(P = P, Q = Q, nrounds = nrounds, ...)
  new("mbtransfer_model", parameters = fit, method = "mbtransfer", hyper = hyper)
}

#' @export
mbtransfer_predict <- function(object, newdata) {
  lags <- time_lags(object@parameters[[1]])
  result <- list()
  subject <- subject_data(newdata) |>
    select(-subject) |>
    as.matrix()
  
  for (i in seq_along(newdata)) {
    result[[i]] <- model_predict_single(
      object@parameters, 
      newdata[[i]],
      lags, 
      subject[i,, drop = FALSE]
    )
  }

  new("ts_inter", series = result, subject_data = subject_data(newdata))
}

model_predict_single <- function(fit, ts_inter, lags, subject = NULL) {
  n_time <- ncol(ts_inter)
  w <- interventions(ts_inter)
  while(ncol(ts_inter) < ncol(w)) {
    ts_inter <- model_predict_step(ts_inter, fit, lags, subject)
  }

  colnames(values(ts_inter)) <- colnames(w)
  ts_inter
}

model_predict_step <- function(ts_inter, fit, lags, subject = NULL) {
  xz <- predictors(ts_inter, lags, subject)
  y_hat <- vector(length = nrow(ts_inter))

  for (j in seq_len(nrow(ts_inter))) {
    y_hat[j] <- predict(fit[[j]], xz)
  }

  values(ts_inter) <- cbind(values(ts_inter), y_hat)
  ts_inter@time <- c(ts_inter@time, max(ts_inter@time) + 1)
  ts_inter
}

#' @export
setClass(
  "mbtransfer_model",
  slots = c(
    parameters = "ANY",
    method = "character",
    hyper = "list"
  )
)

#' @export
setMethod("predict",  c(object = "mbtransfer_model"), mbtransfer_predict)
