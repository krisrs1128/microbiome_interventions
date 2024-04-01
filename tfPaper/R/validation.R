#' @export
train_test_split <- function(ts_inter, p_test = 0.2) {
  N <- length(ts_inter)
  test_ix <- seq_len(N) %in% sample(N, N * p_test)
  list(train = ts_inter[!test_ix], test = ts_inter[test_ix])
}

intervention_start <- function(interventions, offset = 0) {
  perturbed <- colSums(interventions) > 0
  if (any(head(perturbed, offset - 1))) {
    start_ix <- min(which(perturbed))
  } else {
    start_ix <- ncol(interventions) / 2
  }

  start_ix + offset
}

#' @importFrom dplyr bind_rows
#' @export
evaluation <- function(y, y_hat, test_ix) {
  residuals <- list()
  avg_err <- list()

  for (i in seq_along(y)) {
    test_ix[[i]] <- intersect(test_ix[[i]], seq_len(ncol(y_hat[[i]])))
    if (length(test_ix[[i]]) == 0) {
      residuals[[i]] <- rep(NA, nrow(y_hat[[i]]))
    } else {
      residuals[[i]] <- values(y[[i]][, test_ix[[i]]]) - values(y_hat[[i]][, test_ix[[i]]])
    }

    avg_err[[i]] <- data.frame(
      mse = mean(residuals[[i]]^2),
      mae = mean(abs(residuals[[i]])),
      mse_std = mean(residuals[[i]]^2) / var(as.numeric(values(y[[i]][, -test_ix[[i]]]))),
      mae_std = mean(abs(residuals[[i]])) / IQR(as.numeric(values(y[[i]][, -test_ix[[i]]])))
    )
  }

  bind_rows(avg_err, .id = "subject")
}

# #' @importFrom mbtransfer predict
# #' @importFrom mdsine predict
# #' @importFrom fido2 predict
#' @importFrom dplyr bind_rows
#' @importFrom tidyr separate
#' @export
cross_validate <- function(ts_inter, train, K = 5, n_ahead = 5, offset = -1, seeds = NULL) {
  N <- length(ts_inter)

  fits <- list()
  y_hat <- list()
  metrics <- list()

  for (k in seq_len(K)) {
    if (!is.null(seeds)) {
      set.seed(seeds[k])
    }
    splits <- train_test_split(ts_inter, 1 / K)
    fits[[k]] <- train(splits$train)

    stest <- splits$test
    t_starts <- vector(length = length(stest))
    for (i in seq_along(stest)) {
      t_starts[i] <- intervention_start(interventions(stest[[i]]), offset)
      values(stest[[i]]) <- values(stest[[i]])[, seq_len(t_starts[i]), drop = FALSE]

      predict_ix <- seq_len(min(t_starts[i] + n_ahead, ncol(splits$test[[i]])))
      interventions(stest[[i]]) <- interventions(stest[[i]])[, predict_ix, drop = FALSE]
      stest[[i]]@time <- stest[[i]]@time[predict_ix]
    }

    y_hat[[k]] <- predict(fits[[k]], stest)
    for (h in seq_len(n_ahead)) {
      test_ix <- as.list(t_starts + h)
      metrics[[glue("{k}-{h}")]] <- evaluation(splits$test, y_hat[[k]], test_ix)
    }
  }

  metrics <- bind_rows(metrics, .id = "fold_lag") |>
    separate(fold_lag, c("fold", "lag"))
  list(fits = fits, y_hat = y_hat, metrics = metrics)
}
