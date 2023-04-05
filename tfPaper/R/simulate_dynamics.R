
#' @export
step_t <- function(f, g, h) {
  fun <- function(x, w, z) {
    f$fun(x) + g$fun(w) + h$fun(z, w)
  }
  params <- list(P = f$lag, Q = g$lag, PQ = h$lag)
  list(fun = fun, params = params)
}

#' @export
matnorm <- function(N, M, mu = 0, sigma = 1) {
  matrix(rnorm(N * M, mu, sigma), N, M)
}

#' @export
matunif <- function(N, M, lower = -1, upper = 1) {
  matrix(runif(N * M, lower, upper), N, M)
}

#' @export
matunifband <- function(N, M, lower = .5, upper = 1) {
  x <- matrix(runif(N * M, lower, upper), N, M)
  z <- matrix(sample(c(1, -1), N * M, replace = TRUE), N, M)
  x * z
}


#' @importFrom MCMCpack rdirichlet
factorized_step <- function(J, K, D = NULL, ...) {
  if (is.null(D)) {
    D <- J
  }
  
  Theta <- matunif(J, K, ...)
  B <- matunifband(K, D, ...)
  list(B = B, Theta = Theta)
}

#' @export
sparsify <- function(A, sparsity = 0.8) {
  for (i in seq_along(A)) {
    ix <- sample(length(A[[i]]), length(A[[i]]) * sparsity)
    A[[i]][ix] <- 0
  }
  A
}

#' @export
sparsify_rows <- function(B, nonnull_ix = NULL, sparsity_frac = 0.8) {
  if (is.null(nonnull_ix)) {
    nonnull_ix <- sample(nrow(B[[1]]), sparsity_frac * nrow(B[[1]]))
  }
  
  for (i in seq_along(B)) {
    B[[i]][-nonnull_ix, ] <- 0
  }
  B
}

#' @export
linear_sum <- function(coefs) {
  n_lag <- length(coefs)
  n_taxa <- nrow(coefs[[1]])
  
  fun <- function(x) {
    result <- matrix(0, n_taxa, 1)

    L <- ncol(x)
    for (l in seq_len(n_lag)) {
      result <- result + coefs[[l]] %*% x[, L - l + 1]
    }

    result
  }
  
  list(fun = fun, lag = n_lag, coefs = coefs)
}

#' @export
low_rank_step <- function(n_taxa, n_latent, n_perturb = NULL, n_lag = 1, ...) {
  coefs <- list()
  for (l in seq_len(n_lag)) {
    coefs[[l]] <- factorized_step(n_taxa, n_latent, n_perturb, ...)
  }
  
  coefs |>
    map(~ .$Theta %*% .$B)
}

#' @export
normalize_mat <- function(A, lambda = 0.9) {
  for (i in seq_along(A)) {
    A[[i]] <- lambda * A[[i]] / max(svd(A[[i]])$d[1])
  }
  A
}

#' @export
low_rank_step_ <- function(n_taxa, n_latent, n_perturb, n_lag) {
  coefs <- list()
  for (j in seq_len(n_taxa)) {
    coefs[[j]] <- low_rank_step(n_taxa, n_latent, n_perturb, n_lag)
  }
  coefs
}

#' @export
random_step <- function(n_taxa, n_perturb, n_lag, lower = .1, upper = 1) {
  result <- list()
  for (l in seq_len(n_lag)) {
    result[[l]] <- matunifband(n_taxa, n_perturb, lower, upper)
  }
  result
}

#' @export
random_interventions <- function(n_perturb, n_time, n_lag = 3) {
  w <- replicate(n_subject, matrix(0, n_perturb, n_time), simplify = FALSE)
  for (i in seq_along(w)) {
    for (j in seq_len(n_perturb)) {
      start_ix <- sample((n_time / 3) : (n_time - n_lag), 1)
      end_ix <- start_ix + sample((0.5 * n_lag) : (1.5 * n_lag), 1)
      w[[i]][j, start_ix:min(n_time, end_ix)] <- 1
    }
  }
  
  w
}

#' @export
interaction_sum <- function(coefs) {
  n_lag <- length(coefs)
  n_taxa <- nrow(coefs[[1]])
  
  fun <- function(z, w) {
    result <- matrix(0, n_taxa, 1)
    for (l in seq(n_lag)) {
      result <- result + (coefs[[l]] * z) %*% w[, n_lag - l + 1]
    }
    result
  }
  
  list(fun = fun, PQ = n_lag, coefs = coefs)
}

#' @export
gaussian_sampler <- function(sigma = 1) {
  function(theta) {
    if (length(sigma) == 1) {
      sigma <- rep(sigma, length(theta))
    }
  
    rnorm(length(theta), theta, sigma)
  }
}

#' @export
nbinom_sampler <- function(size = 1, baseline = 1) {
  function(log_theta) {
    if (length(size) == 1) {
      size <- rep(size, length(log_theta))
    }
    if (length(baseline) == 1) {
      baseline <- rep(baseline, length(log_theta))
    }

    mu <- baseline * exp(log_theta)
    rnbinom(length(mu), mu  = mu, size = size)
  }
}

#' @export
generate_sample <- function(theta0, w, z, step_generator, sampler) {
  result <- list()
  n_time <- ncol(w)
  n_taxa <- nrow(theta0)
  P <- step_generator$params$P
  Q <- step_generator$params$Q
  
  theta <- cbind(theta0, matrix(0, n_taxa, n_time - ncol(theta0)))
  x <- matrix(0, n_taxa, n_time)
  for (i in seq(max(P, Q) + 1, n_time)) {
    theta[, i] <- step_generator$fun(theta[, (i - P) : (i - 1), drop=F], w[, (i - Q + 1) : i, drop=F], z)
    x[, i] <- sampler(theta[, i])
  }
  
  x
}

#' @export
generate_samples <- function(step_generator, sampler, theta0, w, z) {
  x <- list()
  for (i in seq_len(nrow(z))) {
    x[[i]] <- generate_sample(theta0, w[[i]], z[i, ], step_generator, sampler)
    rownames(x[[i]]) <- str_c("tax", seq_len(nrow(x[[i]])))
  }
  x  
}
