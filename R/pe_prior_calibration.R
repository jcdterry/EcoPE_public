make_known_pe_sources <- function(
    seed = 20260610,
    n = 50000,
    low_cycles = 500,
    low_noise_sd = 0) {
  set.seed(seed)

  low <- sin(seq(0, 2 * pi * low_cycles, length.out = n)) +
    rnorm(n, 0, low_noise_sd)
  high <- rnorm(n)

  list(
    low = low,
    high = high,
    true_pe = c(
      low = normalized_entropy(motif_counts_3(low)),
      high = normalized_entropy(motif_counts_3(high))
    ),
    parameters = c(
      n = n,
      low_cycles = low_cycles,
      low_noise_sd = low_noise_sd
    )
  )
}

make_short_series_count_data <- function(
    source,
    n_series = 10,
    series_length = 50,
    seed = 1,
    random_windows = FALSE) {
  set.seed(seed)

  if (random_windows) {
    starts <- sort(sample(seq_len(length(source) - series_length + 1), n_series))
  } else {
    starts <- seq(1, by = series_length, length.out = n_series)
  }

  rows <- vector("list", n_series)
  for (i in seq_len(n_series)) {
    x <- source[starts[i]:(starts[i] + series_length - 1)]
    counts <- motif_counts_3(x)
    rows[[i]] <- data.frame(
      series_id = i,
      start_index = starts[i],
      end_index = starts[i] + series_length - 1,
      motif_123 = counts[["123"]],
      motif_213 = counts[["213"]],
      motif_312 = counts[["312"]],
      motif_132 = counts[["132"]],
      motif_231 = counts[["231"]],
      motif_321 = counts[["321"]],
      motif_total = sum(counts),
      pe_plugin = normalized_entropy(counts)
    )
  }

  do.call(rbind, rows)
}

make_nsb_baseline_stan_data <- function(count_data, prior_only = 0) {
  y <- as.matrix(count_data[, paste0("motif_", motif_labels_3)])

  list(
    N = nrow(y),
    K = ncol(y),
    y = y,
    motif_total = as.integer(rowSums(y)),
    prior_only = as.integer(prior_only)
  )
}

make_change_gradient_trend_counts <- function(
    seed = 2026,
    n_series = 100,
    n_steps = 50,
    window_size = 25,
    step_size = 5,
    ...) {
  sim <- simulate_change_in_pe_gradient(
    seed = seed,
    n_series = n_series,
    n_steps = n_steps,
    ...
  )

  windows <- make_motif_window_data(
    abundance = sim$abundance,
    window_size = window_size,
    step_size = step_size
  )
  windows$source <- "variance_gradient"

  list(
    windows = windows,
    simulation = sim
  )
}

make_time_varying_known_pe_counts <- function(
    sources,
    n_low = 10,
    n_high = 10,
    series_length = 50,
    seed = 1,
    random_windows = FALSE) {
  low_counts <- make_short_series_count_data(
    source = sources$low,
    n_series = n_low,
    series_length = series_length,
    seed = seed,
    random_windows = random_windows
  )
  high_counts <- make_short_series_count_data(
    source = sources$high,
    n_series = n_high,
    series_length = series_length,
    seed = seed + 1000,
    random_windows = random_windows
  )

  low_counts$source <- "low"
  high_counts$source <- "high"
  out <- rbind(low_counts, high_counts)
  out$period_id <- seq_len(nrow(out))
  out$window_midpoint <- out$period_id
  out$time_scaled <- as.numeric(scale(out$period_id, center = TRUE, scale = TRUE))
  row.names(out) <- NULL
  out
}

make_nsb_time_trend_stan_data <- function(
    count_data,
    prior_only = 0,
    slope_prior_sd = 0.75) {
  y <- as.matrix(count_data[, paste0("motif_", motif_labels_3)])

  list(
    N = nrow(y),
    K = ncol(y),
    y = y,
    time_scaled = count_data$time_scaled,
    motif_total = as.integer(rowSums(y)),
    prior_only = as.integer(prior_only),
    slope_prior_sd = slope_prior_sd
  )
}

summarise_nsb_time_trend_fit <- function(fit, count_data) {
  draws <- as_fit_draw_matrix(
    fit,
    variables = c(
      "pe_start", "pe_end", "pe_change",
      "normalized_pe", "rho_overdisp", "alpha0", "alpha_nsb"
    )
  )

  pe_window <- summarise_draw_matrix(draws, "normalized_pe")
  pe_window <- cbind(
    count_data,
    pe_model_mean = pe_window$mean,
    pe_model_q05 = pe_window$q05,
    pe_model_q50 = pe_window$q50,
    pe_model_q95 = pe_window$q95
  )
  row.names(pe_window) <- NULL

  pe_change <- data.frame(
    variable = c("pe_start", "pe_end", "pe_change"),
    mean = colMeans(draws[, c("pe_start", "pe_end", "pe_change"), drop = FALSE]),
    q05 = apply(draws[, c("pe_start", "pe_end", "pe_change"), drop = FALSE], 2, quantile, probs = 0.05),
    q50 = apply(draws[, c("pe_start", "pe_end", "pe_change"), drop = FALSE], 2, quantile, probs = 0.50),
    q95 = apply(draws[, c("pe_start", "pe_end", "pe_change"), drop = FALSE], 2, quantile, probs = 0.95),
    row.names = NULL
  )

  diagnostics <- data.frame(
    rho_overdisp_mean = mean(draws[, "rho_overdisp"]),
    alpha0_median = unname(quantile(draws[, "alpha0"], 0.50)),
    alpha_nsb_median = unname(quantile(draws[, "alpha_nsb"], 0.50))
  )

  list(
    window_pe = pe_window,
    pe_change = pe_change,
    diagnostics = diagnostics
  )
}

summarise_nsb_event_time_fit <- function(fit, event_data) {
  draws <- as_fit_draw_matrix(
    fit,
    variables = c("pe_start", "pe_end", "pe_change", "pe_time", "alpha_nsb")
  )

  time_lookup <- unique(event_data[, c("time_id", "center_time", "time_scaled")])
  time_lookup <- time_lookup[order(time_lookup$time_id), ]

  pe_time <- summarise_draw_matrix(draws, "pe_time")
  pe_time <- cbind(
    time_lookup,
    pe_model_mean = pe_time$mean,
    pe_model_q05 = pe_time$q05,
    pe_model_q50 = pe_time$q50,
    pe_model_q95 = pe_time$q95
  )
  row.names(pe_time) <- NULL

  pe_change <- data.frame(
    variable = c("pe_start", "pe_end", "pe_change"),
    mean = colMeans(draws[, c("pe_start", "pe_end", "pe_change"), drop = FALSE]),
    q05 = apply(draws[, c("pe_start", "pe_end", "pe_change"), drop = FALSE], 2, quantile, probs = 0.05),
    q50 = apply(draws[, c("pe_start", "pe_end", "pe_change"), drop = FALSE], 2, quantile, probs = 0.50),
    q95 = apply(draws[, c("pe_start", "pe_end", "pe_change"), drop = FALSE], 2, quantile, probs = 0.95),
    row.names = NULL
  )

  diagnostics <- data.frame(
    alpha_nsb_median = unname(quantile(draws[, "alpha_nsb"], 0.50))
  )

  list(
    pe_time = pe_time,
    pe_change = pe_change,
    diagnostics = diagnostics
  )
}

summarise_nsb_timepoint_fit <- function(fit, timepoint_data) {
  draws <- as_fit_draw_matrix(
    fit,
    variables = c("pe_start", "pe_end", "pe_change", "pe_time", "alpha_nsb")
  )

  pe_time <- summarise_draw_matrix(draws, "pe_time")
  pe_time <- cbind(
    timepoint_data,
    pe_model_mean = pe_time$mean,
    pe_model_q05 = pe_time$q05,
    pe_model_q50 = pe_time$q50,
    pe_model_q95 = pe_time$q95
  )
  row.names(pe_time) <- NULL

  pe_change <- data.frame(
    variable = c("pe_start", "pe_end", "pe_change"),
    mean = colMeans(draws[, c("pe_start", "pe_end", "pe_change"), drop = FALSE]),
    q05 = apply(draws[, c("pe_start", "pe_end", "pe_change"), drop = FALSE], 2, quantile, probs = 0.05),
    q50 = apply(draws[, c("pe_start", "pe_end", "pe_change"), drop = FALSE], 2, quantile, probs = 0.50),
    q95 = apply(draws[, c("pe_start", "pe_end", "pe_change"), drop = FALSE], 2, quantile, probs = 0.95),
    row.names = NULL
  )

  diagnostics <- data.frame(
    alpha_nsb_median = unname(quantile(draws[, "alpha_nsb"], 0.50))
  )

  list(
    pe_time = pe_time,
    pe_change = pe_change,
    diagnostics = diagnostics
  )
}

summarise_nsb_event_linear_re_fit <- function(fit, event_data) {
  draws <- as_fit_draw_matrix(
    fit,
    variables = c(
      "pe_population_start", "pe_population_end", "pe_population_change",
      "pe_population_slope_scaled", "pe_population_time", "pe_period_time",
      "beta_time", "alpha_nsb", "sigma_site", "sigma_period"
    )
  )

  time_lookup <- unique(event_data[, c("time_id", "center_time", "time_scaled")])
  time_lookup <- time_lookup[order(time_lookup$time_id), ]

  population_pe <- summarise_draw_matrix(draws, "pe_population_time")
  population_pe <- cbind(
    time_lookup,
    pe_model_mean = population_pe$mean,
    pe_model_q05 = population_pe$q05,
    pe_model_q50 = population_pe$q50,
    pe_model_q95 = population_pe$q95
  )
  row.names(population_pe) <- NULL

  period_pe <- summarise_draw_matrix(draws, "pe_period_time")
  period_pe <- cbind(
    time_lookup,
    pe_model_mean = period_pe$mean,
    pe_model_q05 = period_pe$q05,
    pe_model_q50 = period_pe$q50,
    pe_model_q95 = period_pe$q95
  )
  row.names(period_pe) <- NULL

  change_vars <- c(
    "pe_population_start", "pe_population_end",
    "pe_population_change", "pe_population_slope_scaled"
  )
  pe_change <- data.frame(
    variable = change_vars,
    mean = colMeans(draws[, change_vars, drop = FALSE]),
    q05 = apply(draws[, change_vars, drop = FALSE], 2, quantile, probs = 0.05),
    q50 = apply(draws[, change_vars, drop = FALSE], 2, quantile, probs = 0.50),
    q95 = apply(draws[, change_vars, drop = FALSE], 2, quantile, probs = 0.95),
    row.names = NULL
  )

  beta_time <- summarise_draw_matrix(draws, "beta_time")
  motif_slopes <- data.frame(
    motif_id = beta_time$row_index,
    motif_label = motif_labels_3[beta_time$row_index],
    beta_time_mean = beta_time$mean,
    beta_time_q05 = beta_time$q05,
    beta_time_q50 = beta_time$q50,
    beta_time_q95 = beta_time$q95,
    row.names = NULL
  )

  diagnostics <- data.frame(
    alpha_nsb_median = unname(quantile(draws[, "alpha_nsb"], 0.50)),
    sigma_site_median = median(draws[, grep("^sigma_site\\[", colnames(draws), value = TRUE)]),
    sigma_period_median = median(draws[, grep("^sigma_period\\[", colnames(draws), value = TRUE)])
  )

  list(
    population_pe = population_pe,
    period_pe = period_pe,
    pe_change = pe_change,
    motif_slopes = motif_slopes,
    diagnostics = diagnostics
  )
}

summarise_nsb_event_smooth_walk_re_fit <- function(fit, event_data) {
  draws <- as_fit_draw_matrix(
    fit,
    variables = c(
      "pe_population_start", "pe_population_end", "pe_population_change",
      "pe_population_time", "pe_period_time", "alpha_nsb",
      "sigma_smooth", "sigma_site", "sigma_period"
    )
  )

  time_lookup <- unique(event_data[, c("time_id", "center_time", "time_scaled")])
  time_lookup <- time_lookup[order(time_lookup$time_id), ]

  population_pe <- summarise_draw_matrix(draws, "pe_population_time")
  population_pe <- cbind(
    time_lookup,
    pe_model_mean = population_pe$mean,
    pe_model_q05 = population_pe$q05,
    pe_model_q50 = population_pe$q50,
    pe_model_q95 = population_pe$q95
  )
  row.names(population_pe) <- NULL

  period_pe <- summarise_draw_matrix(draws, "pe_period_time")
  period_pe <- cbind(
    time_lookup,
    pe_model_mean = period_pe$mean,
    pe_model_q05 = period_pe$q05,
    pe_model_q50 = period_pe$q50,
    pe_model_q95 = period_pe$q95
  )
  row.names(period_pe) <- NULL

  pe_change <- data.frame(
    variable = c("pe_population_start", "pe_population_end", "pe_population_change"),
    mean = colMeans(draws[, c("pe_population_start", "pe_population_end", "pe_population_change"), drop = FALSE]),
    q05 = apply(draws[, c("pe_population_start", "pe_population_end", "pe_population_change"), drop = FALSE], 2, quantile, probs = 0.05),
    q50 = apply(draws[, c("pe_population_start", "pe_population_end", "pe_population_change"), drop = FALSE], 2, quantile, probs = 0.50),
    q95 = apply(draws[, c("pe_population_start", "pe_population_end", "pe_population_change"), drop = FALSE], 2, quantile, probs = 0.95),
    row.names = NULL
  )

  diagnostics <- data.frame(
    alpha_nsb_median = unname(quantile(draws[, "alpha_nsb"], 0.50)),
    sigma_smooth_median = median(draws[, grep("^sigma_smooth\\[", colnames(draws), value = TRUE)]),
    sigma_site_median = median(draws[, grep("^sigma_site\\[", colnames(draws), value = TRUE)]),
    sigma_period_median = median(draws[, grep("^sigma_period\\[", colnames(draws), value = TRUE)])
  )

  list(
    population_pe = population_pe,
    period_pe = period_pe,
    pe_change = pe_change,
    diagnostics = diagnostics
  )
}

summarise_nsb_event_smooth_re_fit <- function(...) {
  summarise_nsb_event_smooth_walk_re_fit(...)
}
summarise_nsb_baseline_fit <- function(fit, label, true_pe = NA_real_) {
  draws <- as_fit_draw_matrix(
    fit,
    variables = c("normalized_pe", "rho_overdisp", "alpha0", "alpha_nsb")
  )

  data.frame(
    label = label,
    true_pe = true_pe,
    pe_mean = mean(draws[, "normalized_pe"]),
    pe_q05 = unname(quantile(draws[, "normalized_pe"], 0.05)),
    pe_q50 = unname(quantile(draws[, "normalized_pe"], 0.50)),
    pe_q95 = unname(quantile(draws[, "normalized_pe"], 0.95)),
    rho_overdisp_mean = mean(draws[, "rho_overdisp"]),
    alpha0_median = unname(quantile(draws[, "alpha0"], 0.50)),
    alpha_nsb_median = unname(quantile(draws[, "alpha_nsb"], 0.50)),
    row.names = NULL
  )
}
