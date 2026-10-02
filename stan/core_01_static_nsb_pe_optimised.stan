// Core model 01: static NSB-regularised permutation entropy
//
// Purpose
// -------
// This is the simplest model in the sequence. It estimates one fixed motif
// probability vector from a set of observed ordinal-pattern events and converts
// that posterior motif vector into normalised Shannon permutation entropy (PE).
//
// Intended use
// ------------
// Use this model to demonstrate the small-sample problem: the plug-in PE from a
// short time series can be biased, especially when some motifs are rare or
// unobserved. The model uses the same categorical information as the plug-in
// estimator, but regularises the unknown motif probabilities with an NSB-style
// entropy prior.
//
// Data convention
// ---------------
// The R preprocessing layer should convert a time series into one row per motif
// event:
//
//   y[n] = integer motif label for event n, in 1...K
//
// For 3-motifs, K = 6. For a different motif length, K is the number of possible
// ordinal patterns. The Stan model does not need to know the original time
// series length or motif length separately; it only needs K and the observed
// event labels.

functions {
  // NSB mixture weight for a symmetric Dirichlet concentration alpha.
  //
  // A fixed-alpha symmetric Dirichlet prior does not imply a flat prior over
  // entropy. The NSB idea is to mix over alpha with density proportional to the
  // derivative of the expected entropy with respect to alpha. This makes the
  // induced prior over entropy much closer to uniform than a fixed Dirichlet.
  real nsb_weight_lp(real alpha, int k) {
    real d_expected_h = k * trigamma(k * alpha + 1) - trigamma(alpha + 1);
    return log(abs(d_expected_h));
  }
  
  // Softmax probabilities are unchanged if the same constant is added to every
  // logit. Centering removes this unidentified location. It is used more heavily
  // in models 02 and 03, but included here so the three files share vocabulary.
  vector centered_eta(vector eta) {
    return eta - mean(eta);
  }
  
  // Shannon entropy normalised to [0, 1]. A value of 1 means all motifs are
  // equally likely; lower values mean the motif distribution is more structured.
  real normalised_entropy(vector theta) {
    // This version can handle 0s, necessary for the 'pe_plugin' direct motif comparison.
    int k = num_elements(theta);
    real h = 0;
    for (i in 1:k) {
      if (theta[i] > 0) {
        h -= theta[i] * log(theta[i]);
      }
    }
    return h / log(k);
  }
}


data {
  // Number of motif events after preprocessing. For a single time series of
  // length T and motif length m, this will usually be T - m + 1.
  int<lower=1> N;
  
  // Number of possible motif categories. For 3-motifs, K = 6.
  int<lower=2> K;
  
  // Observed motif category for each event.
  array[N] int<lower=1, upper=K> y;
  
  // If prior_only = 1, the likelihood is skipped. This is useful for checking
  // what the NSB prior alone implies before looking at data.
  int<lower=0, upper=1> prior_only;
}

transformed data {
  // These quantities reproduce the naive plug-in estimator inside Stan so the
  // model output can be compared directly with the Bayesian NSB estimate.
  array[K] int observed_motif_count = rep_array(0, K);
  vector[K] plugin_theta;
  real plugin_pe;
  
  // Count the observed frequencies
  for (n in 1:N) {
    observed_motif_count[y[n]] += 1;
  }
  plugin_theta = to_vector(observed_motif_count) / N; // array-to-vector element-wise division
  plugin_pe = normalised_entropy(plugin_theta);
}


parameters {
  // Unknown fixed motif probability vector.
  simplex[K] theta0;
  
  // Symmetric Dirichlet concentration mixed with the NSB entropy weight.
  // The finite bounds keep the prior proper and avoid pathological numerical
  // values without fixing the amount of smoothing in advance.
  real<lower=1e-8, upper=1e4> alpha_nsb;
}

model {
  // NSB baseline: alpha_nsb is not simply given a generic prior. Instead, the
  // target density is incremented by the entropy-flattening NSB Jacobian, and
  // theta0 is then drawn from the corresponding symmetric Dirichlet.
  target += nsb_weight_lp(alpha_nsb, K);
  theta0 ~ dirichlet(rep_vector(alpha_nsb, K));
  
  if (prior_only == 0) {
    // Event-level categorical likelihood. This is equivalent to a multinomial
    // likelihood on motif counts, but the event format is easier to keep
    // consistent with the trend and multi-site models.
    y ~ categorical(theta0);
  }
}

generated quantities {
  // Main Bayesian PE estimate.
  real pe_nsb = normalised_entropy(theta0);
  
  // Naive plug-in PE from the observed motif frequencies. This is fixed by the
  // data and repeated in every posterior draw for convenience.
  real pe_plugin = plugin_pe;
  
  // Difference between the Bayesian posterior draw and the plug-in estimate.
  // Summarising this quantity shows the direction and size of the correction.
  real pe_nsb_minus_plugin = pe_nsb - pe_plugin;
  
  // Expose the motif counts and plug-in probabilities for table-making.
  vector[K] theta_plugin = plugin_theta;
  array[K] int motif_count = observed_motif_count;
}


