
#' @export
tau_summaries <- function(split) {
  int <- do.call(cbind, map(split, ~ rowSums(.$y1 - .$y0)))
  int_abs <- do.call(cbind, map(split, ~ rowSums(abs(.$y1 - .$y0))))
  list(int = int, int_abs = int_abs)
}

#' @export
temporal_mirror <- function(splits) {
  tau <- map(splits, tau_summaries)
  sgn <- rowMeans(sign(tau[[1]]$int * tau[[2]]$int))
  magnitude <- rowMeans(tau[[1]]$int_abs + tau[[2]]$int_abs)
  sgn * magnitude
}

#' @export
hypothetical <- function(ts, perturb_len = 0, n_lag = 5, ...) {
  fit <- train(ts, ...)
  imagined <- ts
  
  for (i in seq_along(ts)) {
    w <- interventions(ts[[i]])
    w_ <- matrix(0, nrow(w), ncol(w))
    colnames(w_) <- colnames(w)
    rownames(w_) <- rownames(w)
    
    # simulate hypothetical forward from the first perturbation
    perturb_ix <- which(colSums(w) != 0)
    if (length(perturb_ix) == 0 || min(perturb_ix) < n_lag) {
      first_inter <- n_lag
    } else {
      first_inter <- min(perturb_ix)
    }

    if (perturb_len > 0) {
      w_[, seq(first_inter, min(first_inter + perturb_len, ncol(w_)))] <- 1
    }
    
    interventions(imagined[[i]]) <- w_
    v <- values(imagined[[i]])
    values(imagined[[i]]) <- v[, seq_len(first_inter)]
  }
  
  predict(fit, imagined)
}

rearrange_counterfactual <- function(counter) {
  splits <- list()
  for (s in seq_along(counter)) {
    split_s <- list()
    for (i in seq_along(counter[[s]]$y0)) {
      split_s[[i]] <- list(
        y1 = values(counter[[s]]$y1[[i]]),
        y0 = values(counter[[s]]$y0[[i]])
      )
    }
    splits[[s]] <-  split_s
  }
  splits
}

split_estimates_ <- function(ts, split_ix, perturb_len = 5, ...) {
  counterfactual <- list(
    list(y1 = hypothetical(ts[split_ix], perturb_len, ...), y0 = hypothetical(ts[split_ix], ...)),
    list(y1 = hypothetical(ts[-split_ix], perturb_len, ...), y0 = hypothetical(ts[-split_ix], ...))
  )
  
  rearrange_counterfactual(counterfactual)
}

#' @export
split_estimates <- function(ts, n_splits = 1, perturb_len = 5, method = "gbm") {
  splits <- list()
  for (s in seq_len(n_splits)) {
    print(str_c("Evaluating split ", s))
    split_ix <- sample(length(ts), 0.5 * length(ts))
    splits[[s]] <- split_estimates_(ts, split_ix, perturb_len, method = method)
  }
  
  splits
}

multiple_split_mirror <- function(multisplits) {
  ms <- list()
  for (s in seq_along(multisplits)) {
    ms[[s]] <- temporal_mirror(multisplits[[s]])
  }
  ms
}
