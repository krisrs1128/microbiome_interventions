data {
  int<lower=0> P; // lags for factors
  int<lower=0> Q; // lags for interventions
  int<lower=0> K; // number of factors
  int<lower=0> D; // number of interventions
  int<lower=0> N; // number of timepoints
  int<lower=0> V; // number of taxa
  matrix[V, N] x; // one subject's data
  matrix[D, N] interventions; // all interventions
}
parameters {
  array[P] matrix[K, K] A;
  array[Q] matrix[K, D] B;
  matrix[K, V] L;
  matrix[K, N] z;
  real<lower=0> sigma_e;
  real<lower=0> sigma_z;
}

model {
  for (k in 1:K) {
    for (n in 1:max(P, Q)) {
      z[k, k] ~ normal(0, sigma_z);
    }
  }
  
  for (n in (max(P, Q) + 1):N) {
    vector[K] mean_z = rep_vector(0, K);
    for (p in 1:P) {
      mean_z += A[p] * z[, n - p];
    }
    
    for (q in 1:Q) {
      mean_z += B[q] * interventions[, n - q];
    }
    
    for (k in 1:K) {
      z[k, n] ~ normal(mean_z, sigma_z);
    }
    
    vector[V] mu = L' * z[, n];
    for (v in 1:V) {
      x[v, n] ~ normal(mu[v], sigma_e);
    }
  }
}
