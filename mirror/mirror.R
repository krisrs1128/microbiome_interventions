fdp_hat <- function(m) {
  m_sort <- sort(abs(m))
  fdp <- tibble(t = m_sort, fdp = 0)

  for (j in seq_along(m_sort)) {
    fdp$fdp[j] <- sum(m < -m_sort[j]) / max(1, sum(m > m_sort[j]))
  }

  fdp
}

tau_q <- function(fdp, q) {
  if (!any(na.omit(fdp$fdp) < q)) {
    return(NA)
  }

  fdp |>
    filter(fdp < q) |>
    slice_min(t) |>
    pull(t) %>%
    .[[1]]
}

inclusion <- function(s_hat) {
  s_hat <- 1.0 * s_hat
  colMeans(s_hat / rowSums(s_hat))
}

consolidate <- function(s_hat, q) {
  I_hat <- inclusion(s_hat)
  ix <- order(I_hat)
  I_sort <- cumsum(I_hat[ix])
  j_star <- max(which(I_sort <= q))
  I_hat > I_hat[ix[j_star]]
}

selections <- function(m, tau) {
  m > tau
}

multiple_data_splitting <- function(ms, q = 0.1) {
  s_hat <- matrix(FALSE, length(ms), length(ms[[1]]))

  for (k in seq_along(ms)) {
    fdp <- fdp_hat(ms[[k]])
    tau <- tau_q(fdp, q)
    s_hat[k, ] <- selections(ms[[k]], tau)
  }

  list(combined = consolidate(s_hat, q), s_hat = s_hat)
}

z_scores <- function(fit) {
  D <- fit$num.independent.variables
   ranks <- fit$variable.importance |>
     rank(x = _)
   qnorm((ranks - 0.5) / D)
}

ranger_mirror <- function(fmla, data) {
  ix <- sample(nrow(data), 0.5 * nrow(data))
  z1 <- ranger(fmla, data[ix, ], importance = "impurity_corrected", num.trees = 1e4) |>
    z_scores()
  z2 <- ranger(fmla, data[-ix, ], importance = "impurity_corrected", num.trees = 1e4) |>
    z_scores()
  sign(z1 * z2) * (abs(z1) + abs(z2))
}
