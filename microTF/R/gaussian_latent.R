
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

guassian_latent_predict <- function(x_obs, interventions, params) {
  params <- summarize_posterior(params)
  z <- gaussian_infer_latent(x_obs, interventions, params)
  z_future <- gaussian_forecast_latent(z, interventions, params$A, params$B, params$sigma_z)
  x_future <- gaussian_forecast_observed(z_future, params$B, params$L, params$sigma_e)
}

#' @importFrom cmdstanr cmdstan_model
#' @export
gaussian_infer_latent <- function(x_obs, interventions, params, n_draws = 100) {
  data_list <- list(
    x = x_obs,
    interventions = interventions,
    A = params$A,
    B = params$B,
    L = params$L,
    K = params$K,
    sigma_e = params$sigma_e,
    sigma_z = params$sigma_z,
    N = ncol(x_obs),
    V = nrow(x_obs),
    P = length(params$A),
    Q = length(params$B),
    D = nrow(interventions)
  )

  model <- cmdstan_model("microTF/inst/gaussian_infer_latent.stan")
  fit <- model$variational(data_list)
  reshape_forecast(fit$draws("z")[1:n_draws, ])
}

#' @export
gaussian_forecast_latent <- function(z, interventions, A, B, sigma_z, H = 10) {
  forecasts <- list()
  for (i in seq_along(z)) {
    forecasts[[i]] <- gaussian_forecast_latent_(z[[i]], interventions, A, B, sigma_z, H)
  }
  forecasts
}

gaussian_forecast_latent_ <- function(zi, interventions, A, B, sigma_z, H) {
  K <- nrow(zi)
  forecast <- cbind(zi, matrix(0, K, H))
  n_obs <- ncol(zi)
  
  for (h in seq_len(H)) {
    for (p in seq_along(A)) {
      forecast[, n_obs + h] <- forecast[, n_obs + h] + A[[p]] %*% forecast[, n_obs + h - p] 
    }
    for (q in seq_along(B)) {
      forecast[, n_obs + h] <- forecast[, n_obs + h] + B[[q]] %*% interventions[, n_obs + h - q] 
    }
    forecast[, n_obs + h] <- forecast[, n_obs + h] + rnorm(K, 0, sigma_z)
  }

  forecast[, -c(1:n_obs)]
}

gaussian_forecast_observed <- function(z, L, sigma_e) {
  D <- nrow(L)
  H <- ncol(z)
  forecast <- matrix(0, D, H)
  
  for (h in seq_len(H)) {
    forecast[, h] <- L %*% z[, h] + rnorm(D, 0, sigma_e)
  }
  
  forecast
}