
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
  z <- gaussian_infer_latent(x_obs, interventions, params)
  z_future <- gaussian_forecast_latent(z, interventions)
  x_future <- gaussian_forecast_observed(z_future)
}

posterior_means <- function(draws, index_names = c("V", "K")) {
  summarise_draws(draws, "mean") %>%
    mutate(index = str_extract(variable, "[0-9,]+")) %>%
    separate(index, index_names, convert = TRUE) %>%
    select(-variable)
}

pivot_list <- function(x, var1, var2) {
  x %>%
    split(.$P, ~ pivot_wider(names_from = any_of(var1), values_from = mean)) %>%
    map(~ select(., -P)) %>%
    map(~ pivot_wider(., names_from = any_of(var2), values_from = mean)) %>%
    map(~ select(., -var1)) %>%
    map(~ as.matrix(.))
}

#' @importFrom posterior summarise_draws
summarize_posterior <- function(params) {
  L <- posterior_means(params$draws("L"), c("K", "V")) %>%
    pivot_wider(names_from = V, values_from = mean) %>%
    arrange(K) %>%
    column_to_rownames("K")
  A <- posterior_means(params$draws("A"), c("P", "K1", "K2")) %>%
    pivot_list("K1", "K2")
  B <- posterior_means(params$draws("B"), c("P", "K", "D")) %>%
    pivot_list("K", "D")
  sigma_e <- summarise_draws(params$draws("sigma_e"), "mean") %>%
    pull(mean)
  sigma_z <- summarise_draws(params$draws("sigma_z"), "mean") %>%
    pull(mean)
  
  list(L = L, A = A, B = B, sigma_e = sigma_e, sigma_z = sigma_z, K = nrow(A[[1]]))
}

gaussian_infer_latent <- function(x_obs, interventions, params) {
  params <- summarize_posterior(params)
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
  model$variational(data_list)
}