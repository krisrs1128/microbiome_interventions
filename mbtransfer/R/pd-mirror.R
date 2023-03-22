
#' @export
consistency_mirror <- function(effects) {
  sgn <- sign(effects[, 1] * effects[, 2])
  magnitude <- rowMeans(abs(effects))
  sgn * magnitude
}

#' @importFrom dplyr mutate bind_rows row_number
#' @importFrom tibble tibble
#' @export
consistency_mirror_multisplit <- function(effects) {
  ms <- list()
  k <- 1
  for (s in seq_along(effects)) {
    for (lag in seq_len(dim(effects[[s]])[3])) {
      ms[[k]] <- tibble(
        m = consistency_mirror(effects[[s]][,, lag]),
        lag = lag,
        multisplit = s
      ) |>
        mutate(taxon = row_number())
      k <- k + 1
    }
  }

  bind_rows(ms)
}

#' @export
pd_generator <- function(ts, lags) {
  x <- patchify_df(ts, lags[1], lags[2])$x
  
  function(fit, partial_x) {
    pnames <- intersect(colnames(partial_x), colnames(x))
    nx <- nrow(partial_x)
    y_hat <- matrix(nrow = length(fit), ncol = nx)

    for (j in seq_along(fit)) {
      x_ <- x
 
      for (i in seq_len(nx)) {
        x_[, pnames] <- partial_x[i, ]
        y_hat[j, i] <- mean(predict(fit[[j]], x_))
      }
    }
    
    y_hat
  }
}

pd_effects <- function(pd, fit, w0, w1) {
  pd(fit, w1) - pd(fit, w0)
}

pd_splits <- function(ts, w0, w1, n_splits, method = "gbm", ...) {
  effects <- replicate(n_splits, array(dim = c(nrow(ts[[1]]), 2, nrow(w0))), simplify = FALSE)

  for (s in seq_len(n_splits)) {
    print(str_c("Evaluating split ", s))
    split_ix <- sample(length(ts), 0.5 * length(ts))
    fits <- list(
      train(ts[split_ix], method = method, ...)@parameters,
      train(ts[-split_ix], method = method, ...)@parameters
    )

    lags <- time_lags(fits[[1]][[1]])
    pd_fun <- list(
      pd_generator(ts[split_ix], lags),
      pd_generator(ts[-split_ix], lags)
    )
    
    for (i in seq_along(fits)) {
      effects[[s]][,i,] <- pd_effects(pd_fun[[i]], fits[[i]], w0, w1)
    }
  }
  
  effects
}