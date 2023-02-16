
sparse_normal <- function(n, m, mu = 0, sigma = 1, s = .2) {
  x <- matrix(rnorm(n * m, mu, sigma), n, m)
  x[sample(n * m, s * n * m)] <- 0
  x
}

#' @export
gaussian_latent_params <- function(k_factors, n_species, sigma_a, sigma_b, P, Q) {
  A <- map(seq_len(P), ~ sparse_normal(k_factors, k_factors, 0, sigma_a))
  B <- map(seq_len(Q), ~ sparse_normal(k_factors, nrow(interventions), 0, sigma_b))
  L <- matrix(rnorm(n_species * k_factors), n_species, k_factors)
  list(A = A, B = B, L = L)
}

#' @export
gaussian_latent <- function(interventions, params, sigma_z, sigma_e) {
  attach(params)
  k_factors <- ncol(A[[1]])
  n_species <- nrow(L)
  
  z <- matrix(0, nrow = k_factors, ncol = ncol(interventions))
  x <- matrix(0, nrow = n_species, ncol = ncol(interventions))
  for (timepoint in seq_len(ncol(interventions))) {
    if (timepoint > max(length(A) + 1, length(B) + 1)) {
      for (p in seq_along(A)) {
        z[, timepoint] <- z[, timepoint] + A[[p]] %*% z[, timepoint - p]
      }
      for (q in seq_along(B)) {
        z[, timepoint] <- z[, timepoint] + B[[q]] %*% interventions[, timepoint - q, drop = F]
      }
      z[, timepoint] <- z[, timepoint] + rnorm(k_factors, 0, sigma_z)
      
      x[, timepoint] <- L %*% (z[, timepoint, drop = F] + rnorm(k_factors, 0, sigma_e))
    }
  }
  
  params <- list(A = A, B = B, L = L)
  list(z = z, x = x, params = params)
}

#' @export
gaussian_latent_ensemble <- function(interventions, n_subjects = 10, n_species = 250, 
                         k_factors = 3, sigma_a = 0.3, sigma_b = 0.2, 
                         sigma_z = 0.1, sigma_e = 0.2, P = 1, Q = 1) {
  params <- gaussian_latent_params(k_factors, n_species, sigma_a, sigma_b, P, Q)
  ensemble <- list()
  z <- list()

  for (i in seq_len(n_subjects)) {
    samples <- gaussian_latent(interventions, params, sigma_z, sigma_e)
    ensemble[[i]] <- new(
      "ts_inter_single", 
      values = samples$x,
      time = ncol(interventions),
      interventions = interventions
    )
    z[[i]] <- samples$z
  }
  
  list(x = new("ts_inter", series = ensemble), params = params, z = z)
}