// Core model 02: single-site linear trend in NSB-regularised PE
// Built for cmdstan 2.39.0
// NB Certain optimisation features require at least Stan 2.36.
// This is currently ( summer 2026) *not* hosted on the CRAN rstan version

// Purpose
// -------
// This model extends core_01_static_nsb_pe.stan by allowing the motif
// probability vector to change linearly through time. It is intended for a
// single reasonably long time series, where each motif event has a time value.
// It is not build to handle gaps in the data.
//
// ------------------------------
// The model does not fit "PE ~ time" directly. Instead, it fits a linear trend
// in the motif probabilities on the logit scale, then derives PE from those
// probabilities. 
// This is because entropy is a nonlinear summary of a
// full motif distribution: many different motif shifts can produce the same PE.
//
// Data input
// ---------------
// The preprocessing layer should create one row per motif event:
//
//   y[n]           = integer motif label for event n, in 1...K
//   time_scaled[n] = scaled centre time for that motif event
//
// time_scaled should usually be centered near zero and scaled to a sensible
// range. With centered time, theta0 is interpretable as the mid-series motif
// distribution. 
// It is this theta0 that recieves the NSB scaling. 

functions {
  real nsb_weight_lp(real alpha, int k) {
    real d_expected_h = k * trigamma(k * alpha + 1) - trigamma(alpha + 1);
    return log(abs(d_expected_h));
  }
  
  vector centering(vector eta) {
    return eta - mean(eta);
  }
  
  real normalised_entropy(vector theta) {
    int k = num_elements(theta);
    real h = -dot_product(theta, log(theta));  // Using fast vectorised operations
    return h / log(k);
  }
}

data {
  // Same first three fields as model 01: event count, motif count, motif labels.
  int<lower=1> N;
  int<lower=2> K;
  array[N] int<lower=1, upper=K> y;
  
  // Time covariate attached to each observed motif event.
  vector[N] time_scaled;
  
  // Prediction grid for reporting the fitted PE trajectory. For simple use,
  // this can be the sorted unique event centre times. The first and last values
  // define pe_start, pe_end, and pe_change.
  int<lower=1> G;
  vector[G] time_grid_scaled;
  
  int<lower=0, upper=1> prior_only;
  
  // Prior standard deviation for motif-specific time slopes. This should be
  // treated as a prior-sensitivity parameter: it controls how far motif
  // probabilities are allowed to move across the scaled time range.
  real<lower=0> slope_prior_sd;
}


parameters {
  // Baseline motif distribution at time_scaled = 0.
  simplex[K] theta0;
  real<lower=1e-8, upper=1e4> alpha_nsb;
  
  // Beta = Motif-specific linear slopes on the logit scale. The centering
  // the identified quantity because adding the same slope to all motif logits
  // leaves the softmax probabilities unchanged.
  // In STAN this is more efficently done by natively enforcing the sum-to-zero constraint
  sum_to_zero_vector[K] beta_time; 
}


model {
  // Same NSB baseline as model 01.
  target += nsb_weight_lp(alpha_nsb, K);
  
  // Prior on theta, using dynamic alpha
  theta0 ~ dirichlet(rep_vector(alpha_nsb, K));
  
  // Regularsation of linear motif-logit time effect.
  beta_time ~ normal(0, slope_prior_sd);
  
  if (prior_only == 0) {
    // Uses STAN's optimised GLM function that includes the softmax internally
    // https://mc-stan.org/docs/functions-reference/bounded_discrete_distributions.html#categorical-logit-glm
    // Construct a design matrix out of time_scaled, and provide intercept offsets via log(theta0)
    matrix[N, 1] X = to_matrix(time_scaled, N, 1);
    vector[K] alpha = log(theta0);
    matrix[1, K] beta_matrix = to_matrix(beta_time');
    
    y ~ categorical_logit_glm(X, alpha, beta_matrix);
  }
} 

generated quantities {
  // Posterior PE evaluated on the supplied time grid.
  vector[G] pe_time;
  array[G] vector[K] theta_time;
  for (g in 1:G) {
    theta_time[g] = softmax(log(theta0) + beta_time * time_grid_scaled[g]);
    pe_time[g] = normalised_entropy(theta_time[g]);
  }
}


