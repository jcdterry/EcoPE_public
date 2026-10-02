// Core model 03: multi-site linear trend in NSB-regularised PE
// Built for stan 2.39.0
//
// Purpose
// -------
// This model extends core_02_single_site_linear_pe_trend.stan to short,
// replicated time series from multiple sites. It estimates one population-level
// linear trend in motif probabilities while accounting for non-independence
// using motif-specific site effects and motif-specific shared time/period shocks.
//
// Why this structure?
// -------------------
// The data are still individual motif events, not overlapping sliding-window
// counts. Multiple sites provide replication at each motif centre time. The site
// and period effects allow the model to acknowledge that those motif events are
// not all independent: some sites may consistently favour particular motifs,
// and some years/periods may push many sites toward the same motifs.
//
// Population PE
// -------------
// The main estimand is pe_population_time: entropy of the population-level motif
// probabilities after removing realized site and period deviations. pe_period_time
// includes the shared period shock and is useful diagnostically, but it is a
// noisier realised-time quantity rather than the clean population trend.

functions {
  real nsb_weight_lp(real alpha, int k) {
    real d_expected_h = k * trigamma(k * alpha + 1) - trigamma(alpha + 1);
    return log(abs(d_expected_h));
  }
  
  real normalized_entropy(vector theta) {
    int k = num_elements(theta);
    real h = -dot_product(theta, log(theta)); 
    return h / log(k);
  }
  
  matrix double_center(matrix x) {
    
    // Double-centering makes a random-effect matrix easier to interpret.
    // Rows are motifs and columns are sites or periods. After double-centering:
    //   * each motif has average random effect zero across columns,
    //   * each column has average logit shift zero across motifs.
    // The first property keeps the population baseline/trend from being silently
    // absorbed into the random effects. The second removes the softmax location
    // degree of freedom, since adding a common constant to all motif logits within
    // one site or period has no effect on probabilities.
    
    int n_row = rows(x);
    int n_col = cols(x);
    
    // Calculate row averages natively using matrix multiplication
    vector[n_row] row_mean = (x * rep_vector(1.0, n_col)) / n_col;
    
    // Calculate column averages natively using matrix multiplication
    row_vector[n_col] col_mean = (rep_row_vector(1.0, n_row) * x) / n_row;
    
    // Calculate the overall grand mean of all cells
    real grand_mean = mean(x);
    
    // Combine layers using Stan's matrix broadcasting layout
    return x - rep_matrix(row_mean, n_col) - rep_matrix(col_mean, n_row) + grand_mean;
  }
}


data {
  // Same event-level core as models 01 and 02.
  int<lower=1> N;
  int<lower=2> K;
  array[N] int<lower=1, upper=K> y;
  
  // Time covariate for each motif event and grid for reporting the trend.
  vector[N] time_scaled;
  int<lower=1> G;
  vector[G] time_grid_scaled;
  
  // Site and period indexes for each motif event.
  int<lower=1> S;
  array[N] int<lower=1, upper=S> site;
  array[N] int<lower=1, upper=G> time_id;
  
  int<lower=0, upper=1> prior_only;
  
  // Prior scales. These are intentionally data inputs rather than constants so
  // the same compiled model can be used for prior sensitivity checks.
  real<lower=0> slope_prior_sd;
  real<lower=0> sigma_site_scale;
  real<lower=0> sigma_period_scale;
}

parameters {
  // Population baseline and its NSB concentration parameter.
  simplex[K] theta0;
  real<lower=1e-8, upper=1e4> alpha_nsb;
  
  // Population-level linear trend in motif logits.
  vector[K] beta_time_raw;
  
  // Site effects: persistent motif tendencies at particular sites.
  vector<lower=0>[K] sigma_site;
  matrix[K, S] z_site;
  
  // Period effects: shared motif shocks at particular centre times.
  vector<lower=0>[K] sigma_period;
  matrix[K, G] z_period;
}

transformed parameters {
  matrix[K, S] site_effect;
  matrix[K, G] period_effect;
  // Non-centered parameterisation: z_* are standard normal, sigma_* set motif-
  // specific scales. Double-centering then imposes the interpretability
  // constraints described in the function comment above.
  site_effect = double_center(diag_pre_multiply(sigma_site, z_site));
  period_effect = double_center(diag_pre_multiply(sigma_period, z_period));
}

model {
  // Same NSB baseline as models 01 and 02.
  target += nsb_weight_lp(alpha_nsb, K);
  theta0 ~ dirichlet(rep_vector(alpha_nsb, K));
  
  // Population trend prior.
  beta_time_raw ~ normal(0, slope_prior_sd);
  
  // Half-normal scale priors plus standard-normal non-centered effects.
  sigma_site ~ normal(0, sigma_site_scale);
  to_vector(z_site) ~ normal(0, 1);
  sigma_period ~ normal(0, sigma_period_scale);
  to_vector(z_period) ~ normal(0, 1);
  
  if (prior_only == 0) {
    
    //Center the beta_times
    vector[K] beta_time = beta_time_raw - mean(beta_time_raw);
    
    // Compute the log-odds linear predictor vector for each individual observation
    // For every individual observation (N), what is its baseline category prob + its site diff + its year shock?
    
    for (n in 1:N) {  // loop over the observations to handle the heierachical structure smoothly
    vector[K] eta = log(theta0) 
                  + beta_time * time_scaled[n] 
                  + site_effect[, site[n]] 
                  + period_effect[, time_id[n]];
    
    // // Uses STAN's function that includes the softmax internally
    // // https://mc-stan.org/docs/functions-reference/bounded_discrete_distributions.html#categorical-distribution  
    y[n] ~ categorical_logit(eta); 
    }
  }
}

generated quantities {
  vector[G] pe_population_time;
  vector[G] pe_period_time;
  array[G] vector[K] theta_population_time;
  array[G] vector[K] theta_period_time;
  
  // Posterior regression parameters
  real pe_slope;
  real pe_intercept;
  real mean_t;
  real mean_H;
  real num;
  real den;
  
  // Center the slope vector locally here too for the prediction grid math
  vector[K] beta_time = beta_time_raw - mean(beta_time_raw);
  
  for (g in 1:G) {
    vector[K] eta_population = log(theta0) + beta_time * time_grid_scaled[g];
    vector[K] eta_period = eta_population + period_effect[, g];
    
    theta_population_time[g] = softmax(eta_population); 
    theta_period_time[g] = softmax(eta_period);
    
    pe_population_time[g] = normalized_entropy(theta_population_time[g]);
    pe_period_time[g] = normalized_entropy(theta_period_time[g]);
  }
  
  // Fit the linear model to this draw's entropy vector using closed-form OLS
  mean_t = mean(time_grid_scaled);
  mean_H = mean(pe_population_time);
  num = 0;
  den = 0;
  for (g in 1:G) {
    real t_dev = time_grid_scaled[g] - mean_t;
    num += t_dev * (pe_population_time[g] - mean_H);
    den += t_dev * t_dev;
  }
  // Intercept and slope for this draw
  pe_slope = num / den;
  pe_intercept = mean_H - pe_slope * mean_t;
  
}
