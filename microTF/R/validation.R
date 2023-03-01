
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
evaluation <- function(y, y_hat, test_ix, window = 10) {
  residuals <- list()
  avg_err <- list()
  
  for (i in seq_along(y)) {
    residuals[[i]]<- values(y[[i]][, test_ix]) - values(y_hat[[i]][, test_ix])
    avg_err[[i]] <- data.frame(
      mse = mean(residuals[[i]] ^ 2),
      mae = mean(abs(residuals[[i]]))
    )
  }
  
  bind_rows(avg_err, .id = "subject")
}

#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @export
cross_validate <- function(ts_inter, train, K = 5, n_ahead = 5, t_start = 5) {
  N <- length(ts_inter)

  fits <- list()
  y_hat <- list()
  metrics <- list()
  
  for (k in seq_len(K)) {
    splits <- train_test_split(ts_inter, 1 / K)
    fits[[k]] <- train(splits$train)

    stest <- splits$test
    for (i in seq_along(splits$test)) {
      values(stest[[i]]) <- values(stest[[i]])[, seq_len(t_start), drop = FALSE]
      interventions(stest[[i]]) <- interventions(stest[[i]])[, seq_len(t_start + n_ahead), drop = FALSE]
    }
    
    y_hat[[k]] <- predict(fits[[k]], stest)
    test_ix <- seq(t_start, t_start + n_ahead)
    metrics[[k]] <- evaluation(splits$test, y_hat[[k]], test_ix, k)
  }
  
  metrics <- bind_rows(metrics, .id = "fold")
  list(fits = fits, y_hat = y_hat, metrics = metrics)
}