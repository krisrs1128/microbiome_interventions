
#' @export
step_t <- function(f, g, h) {
  fun <- function(x, w, z) {
    f$fun(x) + g$fun(w) + h$fun(x, w)
  }
  params <- list(P = f$lag, Q = g$lag, PQ = h$lag)
  list(fun = fun, params = params)
}

#' @export
matnorm <- function(N, M, mu = 0, sigma = 1) {
  matrix(rnorm(N * M, mu, sigma), N, M)
}

#' @importFrom MCMCpack rdirichlet
factorized_step <- function(J, K, D = NULL, ...) {
  if (is.null(D)) {
    D <- J
  }
  
  Theta <- matnorm(J, K, ...)
  B <- matnorm(K, D, ...)
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
  
  coefs %>%
    map(~ .$Theta %*% .$B)
}

#' @export
normalize <- function(A, lambda = 0.9) {
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
interaction_sum <- function(coefs) {
  n_taxa <- length(coefs)
  n_lag <- length(coefs[[1]])
  
  fun <- function(x, w) {
    result <- matrix(0, n_taxa, 1)
 
    for (l in seq(n_lag)) {
      for (j in seq_len(n_taxa)) {
        result[j, 1] <- result[j, 1] + t(x[, n_lag - l + 1]) %*% coefs[[j]][[l]] %*% w[, n_lag - l + 1]
      }
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