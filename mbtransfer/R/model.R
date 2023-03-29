
#' @importFrom xgboost xgboost
#' @export
mbtransfer <- function(ts_inter, P = 1, Q = 1, verbose = 0, nrounds = 25, objective = "reg:pseudohubererror", ...) {
  train_data <- patchify_df(ts_inter, P, Q)
  fit <- list()
  for (j in seq_along(train_data$y)) {
    fit[[j]] <- xgboost(data = train_data$x, label = train_data$y[[j]], nrounds = nrounds, verbose = verbose, objective = objective, ...)
  }
  
  fit
}

#' @export
split_future <- function(series_i) {
  t1 <- ncol(series_i)
  t2 <- ncol(interventions(series_i))
  pre <- series_i[, seq_len(t1), drop = FALSE]
  interventions <- interventions(series_i)[, seq(t1 + 1, t2), drop = FALSE]
  list(pre = pre, interventions = interventions)
}

#' @export
mbtransfer_predict <- function(fit, ts_inter, new_interventions) {
  lags <- time_lags(fit[[1]])
  if (!is.list(new_interventions)) {
    new_interventions <- replicate(length(ts_inter), new_interventions, simplify = FALSE)
  }

  result <- list()
  subject <- subject_data(ts_inter) |>
    select(-subject) |>
    as.matrix()
  
  for (i in seq_along(ts_inter)) {
    result[[i]] <- model_predict_single(
      fit, 
      ts_inter[[i]],
      new_interventions[[i]], 
      lags, 
      subject[i,, drop = FALSE]
    )
  }

  new("ts_inter", series = result)
}

model_predict_single <- function(fit, ts_inter, new_interventions, lags, subject = NULL) {
  for (h in seq_len(ncol(new_interventions))) {
    ts_inter <- model_predict_step(ts_inter, fit, new_interventions[, h, drop = FALSE], lags, subject)
  }
  ts_inter
}

model_predict_step <- function(ts_inter, fit, z_next, lags, subject = NULL) {
  xz <- predictors(ts_inter, z_next, lags, subject)
  y_hat <- vector(length = nrow(ts_inter))
  for (j in seq_len(nrow(ts_inter))) {
    y_hat[j] <- predict(fit[[j]], xz)
  }

  values(ts_inter) <- cbind(values(ts_inter), y_hat)
  ts_inter@time <- c(ts_inter@time, max(ts_inter@time) + 1)
  interventions(ts_inter) <- cbind(interventions(ts_inter), z_next)
  ts_inter
}