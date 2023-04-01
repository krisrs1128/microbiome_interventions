
#' ts is a single element of a ts_inter class
#' @export
replace_inter_ <- function(ts, new_inter, start_ix = NULL) {
  if (is.null(start_ix)) {
    start_ix <- ncol(ts) - ncol(new_inter) + 1
  }

  inter <- interventions(ts)[, seq_len(start_ix), drop = FALSE]
  interventions(ts) <- cbind(inter, new_inter)
  ts
}

#' @export
replace_inter <- function(ts, new_inter, start_ix = NULL) {
  if (length(start_ix) == 1) {
    start_ix <- rep(start_ix, length(ts))
  }
  
  for (i in seq_along(ts)) {
    ts[[i]] <- replace_inter_(ts[[i]], new_inter, start_ix)
  }
  
  ts
}

#' @export
replace_subject <- function(ts, new_subject) {
  subject <- subject_data(ts)
  shared_cols <- intersect(colnames(subject), new_subject)
  subject[, shared_cols] <- new_subject
  subject_data(ts) <- subject
  ts
}

sample_ts <- function(ts, n, patch_len = 5) {
  # randomly subset series
  weights <- map_dbl(ts, ncol)
  weights <- weights / sum(weights)
  ix <- sample(seq_along(ts), n, replace = TRUE, prob = weights)
  ts_star <- ts[ix]
  
  # randomly subset windows
  for (i in seq_along(ts_star)) {
    start_ix <- sample(seq_len(ncol(ts_star[[i]]) - patch_len + 1), 1)
    ts_star[[i]] <- ts_star[[i]][, seq(start_ix, start_ix + patch_len - 1)]
  }
  
  ts_star
}