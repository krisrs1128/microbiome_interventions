
dlm_params <- function(k_factors = 2, n_species = 20, sigma_a = 0.2, P = 1, Q = 1) {
  A <- map(seq_len(P), ~ matrix(rnorm(k_factors ^ 2, 0, sigma_a), k_factors, k_factors))
  B <- map(seq_len(Q), ~ matrix(rnorm(k_factors * nrow(interventions)), k_factors, nrow(interventions)))
  L <- matrix(rnorm(n_species * k_factors), n_species, k_factors)
  list(A = A, B = B, L = L)
}

dlm <- function(interventions, params, sigma_e = 0.2) {
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
      
      x[, timepoint] <- L %*% (z[, timepoint, drop = F] + rnorm(k_factors, 0, sigma_e))
    }
  }
  
  params <- list(A = A, B = B, L = L)
  list(z = z, x = x, params = params)
}

dlm_ensemble <- function(interventions, n_subjects = 10, n_species = 20, k_factors = 2, sigma_a = 0.1) {
  params <- dlm_params(k_factors, n_species, sigma_a)
  ensemble <- list()
  z <- list()

  for (i in seq_len(n_subjects)) {
    samples <- dlm(interventions, params)
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