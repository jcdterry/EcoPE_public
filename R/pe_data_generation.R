motif_labels_3 <- c("123", "213", "312", "132", "231", "321")

motif_counts_3 <- function(x) {
    counts <- statcomp::ordinal_pattern_distribution(x, 3)
    names(counts) <- motif_labels_3
    return(counts)
}

motif_time_series_3 <- function(x) {
    motifs <- statcomp::ordinal_pattern_time_series(x, 3)[seq_len(length(x) - 2)] + 1
    return(as.integer(motifs))
}

normalized_entropy <- function(counts) {
  p <- counts / sum(counts)
  p <- p[p > 0]
  -sum(p * log(p)) / log(length(counts))
}

simulate_feedback_sites <- function(
    seed = 20260609,
    n_site = 100,
    n_year = 50,
    burn_in = 100,             # Added 5-step stabilization period
    delayed_feedback = -1.1,
    r_max = 2.6, 
    K = 10,                  # Carrying capacity (populations oscillate around 10, not 1)
    imm = 0.05,              # Minor immigration floor to prevent true extinction
    sd_start = 0.1,
    sd_end = 0.2,
    rho_spatial = 0.25,
    log_equilibrium = 0) {
  
  set.seed(seed)
  
  # Calculate total simulated timeline length including burn-in steps
  total_steps <- n_year + burn_in
  
  # Initialise tracking matrix 
  abundance <- matrix(K, nrow = total_steps, ncol = n_site)
  abundance[1, ] <- runif(n_site, K*0.4, K*0.6)
  abundance[2, ] <- runif(n_site, K*0.4, K*0.6)
  
  # Pre-calculate standard normal regional shared climate shocks across the timeline
  regional_z <- rnorm(total_steps, 0, 1)
  
  for (t in 3:total_steps) {
    if (t <= (burn_in + 2)) {
      sd_t <- sd_start   # During burn-in, noise variance remains flat
    } else {
      progress <- (t - (burn_in + 2)) / (total_steps - (burn_in + 2))
      sd_t <- sd_start + (sd_end - sd_start) * progress
    }
    
    var_total_t <- sd_t^2
    
    # Partition variance
    sd_shared_t <- sqrt(var_total_t * rho_spatial)
    sd_local_t  <- sqrt(var_total_t * (1 - rho_spatial))
    shared_shock <- regional_z[t] * sd_shared_t
    local_shock  <- rnorm(n_site, 0, sd_local_t)
   
   # RICKER MODEL:
    density_modifier <- (1.0 - (abundance[t - 1, ] / K)) + (delayed_feedback * (abundance[t - 2, ] / K))
    growth_exponent  <- (r_max * density_modifier) + shared_shock + local_shock
    abundance[t, ] <- imm + (abundance[t - 1, ] * exp(growth_exponent))
     }
   
  # Strip burn-in
  post_burn_abundance <- abundance[(burn_in + 1):total_steps, , drop = FALSE]
  colnames(post_burn_abundance) <- paste0("site_", seq_len(n_site))
  
  return(list( abundance_matrix=post_burn_abundance ))
}


make_motif_window_data <- function(
    abundance,
    window_size = 25,
    step_size = 5,
    site_ids = seq_len(ncol(abundance))) {
  n_year <- nrow(abundance)
  starts <- seq(1, n_year - window_size + 1, by = step_size)
  rows <- vector("list", length(starts) * length(site_ids))
  idx <- 1

  for (site in site_ids) {
    for (start in starts) {
      end <- start + window_size - 1
      counts <- motif_counts_3(abundance[start:end, site])
      rows[[idx]] <- data.frame(
        site_id = site,
        period_id = match(start, starts),
        window_start = start,
        window_end = end,
        window_midpoint = (start + end) / 2,
        motif_123 = counts[["123"]],
        motif_213 = counts[["213"]],
        motif_312 = counts[["312"]],
        motif_132 = counts[["132"]],
        motif_231 = counts[["231"]],
        motif_321 = counts[["321"]],
        motif_total = sum(counts),
        pe_plugin = normalized_entropy(counts)
      )
      idx <- idx + 1
    }
  }

  out <- do.call(rbind, rows)
  out$time_scaled <- as.numeric(scale(out$window_midpoint, center = TRUE, scale = TRUE))
  out
}

make_motif_event_data <- function(
    abundance,
    site_ids = seq_len(ncol(abundance))) {
  n_year <- nrow(abundance)
  rows <- vector("list", length(site_ids))

  for (idx in seq_along(site_ids)) {
    site <- site_ids[idx]
    motifs <- motif_time_series_3(abundance[, site])
    rows[[idx]] <- data.frame(
      site_id = site,
      center_time = seq(2, n_year - 1),
      motif_id = motifs,
      motif_label = motif_labels_3[motifs]
    )
  }

  out <- do.call(rbind, rows)
  out$time_id <- as.integer(factor(out$center_time, levels = sort(unique(out$center_time))))
  out$time_scaled <- as.numeric(scale(out$center_time, center = TRUE, scale = TRUE))
  row.names(out) <- NULL
  out
}

make_single_site_stan_data <- function(window_data) {
  y <- as.matrix(window_data[, paste0("motif_", motif_labels_3)])

  list(
    N = nrow(y),
    K = ncol(y),
    y = y,
    time_scaled = window_data$time_scaled,
    motif_total = as.integer(rowSums(y))
  )
}

make_event_time_trend_stan_data <- function(event_data, prior_only = 0, slope_prior_sd = 0.75) {
  time_lookup <- unique(event_data[, c("time_id", "center_time", "time_scaled")])
  time_lookup <- time_lookup[order(time_lookup$time_id), ]

  list(
    N = nrow(event_data),
    K = length(motif_labels_3),
    y = as.integer(event_data$motif_id),
    time_scaled = event_data$time_scaled,
    G = nrow(time_lookup),
    time_grid_scaled = time_lookup$time_scaled,
    prior_only = as.integer(prior_only),
    slope_prior_sd = slope_prior_sd
  )
}

make_event_linear_re_stan_data <- function(
    event_data,
    prior_only = 0,
    slope_prior_sd = 0.35,
    sigma_site_scale = 0.35,
    sigma_period_scale = 0.15) {
  site_index <- as.integer(factor(event_data$site_id))
  time_index <- as.integer(factor(event_data$time_id))
  time_lookup <- unique(event_data[, c("time_id", "center_time", "time_scaled")])
  time_lookup <- time_lookup[order(time_lookup$time_id), ]

  list(
    N = nrow(event_data),
    K = length(motif_labels_3),
    S = max(site_index),
    G = max(time_index),
    y = as.integer(event_data$motif_id),
    site = site_index,
    time_id = time_index,
    time_scaled = event_data$time_scaled,
    time_grid_scaled = time_lookup$time_scaled,
    prior_only = as.integer(prior_only),
    slope_prior_sd = slope_prior_sd,
    sigma_site_scale = sigma_site_scale,
    sigma_period_scale = sigma_period_scale
  )
}

make_event_smooth_walk_re_stan_data <- function(
    event_data,
    prior_only = 0,
    sigma_smooth_scale = 0.20,
    sigma_site_scale = 0.35,
    sigma_period_scale = 0.15) {
  site_index <- as.integer(factor(event_data$site_id))
  time_index <- as.integer(factor(event_data$time_id))

  list(
    N = nrow(event_data),
    K = length(motif_labels_3),
    S = max(site_index),
    G = max(time_index),
    y = as.integer(event_data$motif_id),
    site = site_index,
    time_id = time_index,
    prior_only = as.integer(prior_only),
    sigma_smooth_scale = sigma_smooth_scale,
    sigma_site_scale = sigma_site_scale,
    sigma_period_scale = sigma_period_scale
  )
}

make_event_smooth_re_stan_data <- function(...) {
  make_event_smooth_walk_re_stan_data(...)
}

make_timepoint_count_data <- function(event_data) {
  labels <- motif_labels_3
  time_lookup <- unique(event_data[, c("time_id", "center_time", "time_scaled")])
  time_lookup <- time_lookup[order(time_lookup$time_id), ]

  rows <- vector("list", nrow(time_lookup))
  for (i in seq_len(nrow(time_lookup))) {
    events_i <- event_data[event_data$time_id == time_lookup$time_id[i], ]
    counts <- tabulate(events_i$motif_id, nbins = length(labels))
    names(counts) <- labels
    rows[[i]] <- data.frame(
      time_id = time_lookup$time_id[i],
      center_time = time_lookup$center_time[i],
      time_scaled = time_lookup$time_scaled[i],
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

make_timepoint_multinomial_stan_data <- function(
    timepoint_data,
    prior_only = 0,
    slope_prior_sd = 0.75) {
  y <- as.matrix(timepoint_data[, paste0("motif_", motif_labels_3)])

  list(
    G = nrow(y),
    K = ncol(y),
    y = y,
    time_scaled = timepoint_data$time_scaled,
    prior_only = as.integer(prior_only),
    slope_prior_sd = slope_prior_sd
  )
}

make_multisite_stan_data <- function(window_data) {
  y <- as.matrix(window_data[, paste0("motif_", motif_labels_3)])
  site_index <- as.integer(factor(window_data$site_id))
  period_index <- as.integer(factor(window_data$period_id))

  list(
    N = nrow(y),
    K = ncol(y),
    S = max(site_index),
    P = max(period_index),
    y = y,
    site = site_index,
    period = period_index,
    time_scaled = window_data$time_scaled,
    motif_total = as.integer(rowSums(y))
  )
}

as_abundance_matrix <- function(simulated_data) {
  if (is.list(simulated_data) && "abundance_matrix" %in% names(simulated_data)) {
    abundance <- simulated_data$abundance_matrix
  } else {
    abundance <- simulated_data
  }

  abundance <- as.matrix(abundance)

  if (!is.numeric(abundance)) {
    stop("Abundance data must be numeric.", call. = FALSE)
  }

  if (nrow(abundance) < 3) {
    stop("At least three time points are needed to make 3-point motifs.", call. = FALSE)
  }

  if (is.null(colnames(abundance))) {
    colnames(abundance) <- paste0("site_", seq_len(ncol(abundance)))
  }

  abundance
}

mean_pairwise_site_correlation <- function(simulated_data) {
  abundance <- as_abundance_matrix(simulated_data)

  if (ncol(abundance) < 2) {
    stop("At least two sites are needed for pairwise correlations.", call. = FALSE)
  }

  cor_matrix <- cor(abundance, use = "pairwise.complete.obs")
  mean(cor_matrix[upper.tri(cor_matrix)], na.rm = TRUE)
}

