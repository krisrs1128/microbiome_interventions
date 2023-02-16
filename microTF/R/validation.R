
#' @export
train_test_split <- function(ts_inter, p_test = 0.2) {
  N <- length(ts_inter)
  test_ix <- seq_len(N) %in% sample(N, N * p_test)
  list(train = ts_inter[!test_ix], test = ts_inter[test_ix])
}

#' @export
post_intervention <- function(interventions, window) {
  result <- rep(FALSE, length(interventions))
  for (i in seq_along(result)) {
    if (any(interventions[max(i - window, 1):i] == 1)) {
      result[i] <- TRUE
    }
  }
  
  result
}

#' @importFrom dplyr bind_rows
#' @export
evaluation <- function(y, y_hat, window = 10) {
  residuals <- list()
  avg_err <- list()
  
  for (i in seq_along(y)) {
    residuals[[i]]<- y[[i]]@values - y_hat[[i]]@values
    time_window <- post_intervention(y[[i]]@interventions, window)
    avg_err[[i]] <- data.frame(
      mse = mean(residuals[[i]][, time_window] ^ 2),
      mae = mean(abs(residuals[[i]][, time_window]))
    )
  }
  
  bind_rows(avg_err, .id = "subject")
}

#' @importFrom dplyr bind_rows
#' @export
cross_validate <- function(ts_inter, train, n_ahead = 1, K = 5, n_reps = 5) {
  N <- length(ts_inter)

  fits <- list()
  y_hat <- list()
  metrics <- list()
  
  for (k in seq_len(K)) {
    splits <- train_test_split(ts_inter, 1 / K)
    fits[[k]] <- train(splits$train)
    y_hat[[k]] <- predict(fits[[k]], splits$test, n_ahead)
    metrics[[k]] <- evaluation(splits$test, y_hat[[k]], k)
  }
  
  metrics <- bind_rows(metrics, .id = "fold")
  list(fits = fits, y_hat = y_hat, metrics = metrics)
}