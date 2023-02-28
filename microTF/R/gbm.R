
gbm_predict <- function(object, newdata) {
  fit <- object@parameters
  series <- list()
  new_interventions <- list()
  
  for (i in seq_along(newdata)) {
    t1 <- ncol(newdata[[i]])   
    t2 <- ncol(interventions(newdata[[i]]))
    
    series[[i]] <- ts_inter[[i]][, seq_len(t1)]
    new_interventions[[i]] <- interventions(newdata[[i]])[, seq(t1 + 1, t2)]
  }
  
  ts_inter <- new("ts_inter", series = series)
  gbm_predict_(fit, ts_inter, new_interventions)
}

#' @export
gbm_predict_ <- function(fit, ts_inter, new_interventions) {
  lags <- time_lags(fit[[1]])  
  if (!is.list(new_interventions)) {
    new_interventions <- replicate(length(ts_inter), new_interventions, simplify = FALSE)
  }
  
  result <- list()
  for (i in seq_along(ts_inter)) {
    result[[i]] <- gbm_predict_single(fit, ts_inter[[i]], new_interventions[[i]], lags)
  }
  
  new("ts_inter", series = result)
}

#' @export
gbm_predict_single <- function(fit, ts_inter, new_interventions, lags) {
  for (h in seq_len(ncol(new_interventions))) {
    ts_inter <- gbm_predict_step(ts_inter, fit, new_interventions[, h, drop = FALSE], lags)
  }
  ts_inter
}

#' @export
gbm_predict_step <- function(ts_inter, fit, z_next, lags) {
  xz <- predictors(ts_inter, z_next, lags)
  y_hat <- vector(length = nrow(ts_inter))
  for (j in seq_len(nrow(ts_inter))) {
    y_hat[j] <- predict(fit[[j]], xz)
  }

  values(ts_inter) <- cbind(values(ts_inter), y_hat)
  ts_inter@time <- c(ts_inter@time, max(ts_inter@time) + 1)
  interventions(ts_inter) <- cbind(interventions(ts_inter), z_next)
  ts_inter
}