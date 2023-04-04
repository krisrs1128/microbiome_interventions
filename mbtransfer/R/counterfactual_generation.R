
###############################################################################
# Helpers to generate counterfactual ts_inter objects
###############################################################################

to_vector <- function(x) {
  if (length(x) == 1) {
    x <- c(x)
  }
  x
}

#' Hypothetical Step Interventions
#' @example
#' steps(c("P1" = TRUE), 1:3, 2:3, 4)
#' @export
steps <- function(p_states, starts = 1, lengths = 1:3, L = 3, w_star = c(0, 1)) {
  w0 <- matrix(0, nrow = length(p_states), ncol = L)
  rownames(w0) <- names(p_states)
  colnames(w0) <- glue("T+{seq_len(ncol(w0))}")
  starts <- to_vector(starts)
  lengths <- to_vector(lengths)
  active_p <- names(p_states[p_states])

  result <- list()
  k <- 1
  for (i in seq_along(w_star)) {
    for (s in seq_along(starts)) {
      for (l in seq_along(lengths)) {
        wi <- w0
        for (p in seq_along(active_p)) {
          wi[p, seq(s, min(s + l, ncol(wi)))] <- w_star[i]
        }
        result[[k]] <- wi
        k <- k + 1
      }
    }
  }
  
  unique(result)
}

#' Hypothetical Pulse Interventions
#' @example
#' pulses(c("P1" = TRUE), 1, 4)
#' pulses(c("P1" = TRUE), 1:3, 4)
#" pulses(c("P1" = TRUE), 1:3, 4, seq(0, 1, .2))
#' @export
pulses <- function(p_states, lags = 1, L = 3, w_star = c(0, 1)) {
  w0 <- matrix(0, nrow = nrow(inter), ncol = L)
  rownames(w0) <- rownames(inter)
  colnames(w0) <- glue("T+{seq_len(ncol(w0))}")
  lengths <- to_vector(lengths)
  active_p <- names(p_states[p_states])
  
  result <- list()
  k <- 1
  for (i in seq_along(w_star)) {
    for (l in seq_along(lags)) {
      wi <- w0
      for (p in seq_along(active_p)) {
        wi[p, L - l + 1] <- w_star[i]
      }
      result[[k]] <- wi
      k <- k + 1
    }
  }
  
  result
}

#' Generate Counterfactual versions of a ts_inter object
#' @export
counterfactual_ts <- function(ts, w0, w1) {
  ts0 <- replace_inter(ts, w0)
  ts1 <- replace_inter(ts, w1)
  list(ts0 = ts0, ts1 = ts1)
}
