select_draw_variables <- function(draws, variables) {
  if (is.null(variables)) {
    return(draws)
  }

  selected <- unique(unlist(lapply(variables, function(variable) {
    exact <- intersect(colnames(draws), variable)
    indexed <- grep(paste0("^", variable, "\\["), colnames(draws), value = TRUE)
    c(exact, indexed)
  })))

  draws[, selected, drop = FALSE]
}

load_fit_object <- function(path) {
  if (!is.character(path) || length(path) != 1 || !file.exists(path)) {
    stop("path must be a single existing file")
  }

  rds_obj <- tryCatch(readRDS(path), error = function(e) e)
  if (!inherits(rds_obj, "error")) {
    return(rds_obj)
  }

  env <- new.env(parent = emptyenv())
  loaded <- tryCatch(load(path, envir = env), error = function(e) e)
  if (inherits(loaded, "error")) {
    stop("Could not read fit path as RDS or RData: ", path)
  }
  if (length(loaded) != 1) {
    stop("Expected exactly one object in saved fit file, found: ", paste(loaded, collapse = ", "))
  }

  env[[loaded]]
}

as_fit_draw_matrix <- function(fit, variables = NULL) {
  if (is.character(fit)) {
    fit <- load_fit_object(fit)
  }

  if (inherits(fit, "CmdStanMCMC")) {
    draws <- fit$draws(variables = variables, format = "matrix")
    return(as.matrix(draws))
  }

  if (inherits(fit, "draws")) {
    draws <- posterior::as_draws_matrix(fit)
    return(as.matrix(select_draw_variables(draws, variables)))
  }

  if (is.data.frame(fit) || is.matrix(fit)) {
    draws <- as.matrix(fit)
    return(select_draw_variables(draws, variables))
  }

  stop("fit must be a CmdStanMCMC object, posterior draws object, matrix, or data frame")
}

summarise_draw_matrix <- function(draws, prefix) {
  vars <- grep(paste0("^", prefix, "\\["), colnames(draws), value = TRUE)
  if (!length(vars)) {
    stop("No variables found with prefix: ", prefix)
  }

  index <- as.integer(sub(paste0("^", prefix, "\\[([0-9]+)\\]$"), "\\1", vars))
  vars <- vars[order(index)]
  index <- sort(index)

  data.frame(
    row_index = index,
    mean = colMeans(draws[, vars, drop = FALSE]),
    q05 = apply(draws[, vars, drop = FALSE], 2, quantile, probs = 0.05),
    q50 = apply(draws[, vars, drop = FALSE], 2, quantile, probs = 0.50),
    q95 = apply(draws[, vars, drop = FALSE], 2, quantile, probs = 0.95),
    row.names = NULL
  )
}

softmax_vec <- function(x) {
  z <- x - max(x)
  exp(z) / sum(exp(z))
}

entropy_from_prob <- function(p) {
  p <- p[p > 0]
  -sum(p * log(p)) / log(length(p))
}

extract_parameter_matrix <- function(draws, prefix, k) {
  vars <- paste0(prefix, "[", seq_len(k), "]")
  missing <- setdiff(vars, colnames(draws))
  if (length(missing)) {
    stop("Missing parameter columns: ", paste(missing, collapse = ", "))
  }
  draws[, vars, drop = FALSE]
}

summarise_population_pe_time <- function(fit, window_data, k = 6) {
  draws <- as_fit_draw_matrix(
    fit,
    variables = c(
      paste0("intercept[", seq_len(k), "]"),
      paste0("beta_time[", seq_len(k), "]")
    )
  )

  intercept <- extract_parameter_matrix(draws, "intercept", k)
  beta_time <- extract_parameter_matrix(draws, "beta_time", k)

  time_lookup <- unique(window_data[, c("period_id", "window_start", "window_end",
                                        "window_midpoint", "time_scaled")])
  time_lookup <- time_lookup[order(time_lookup$window_midpoint), ]

  summaries <- vector("list", nrow(time_lookup))
  for (i in seq_len(nrow(time_lookup))) {
    t_scaled <- time_lookup$time_scaled[i]
    pe <- numeric(nrow(draws))

    for (d in seq_len(nrow(draws))) {
      eta <- intercept[d, ] + beta_time[d, ] * t_scaled
      eta <- eta - mean(eta)
      pe[d] <- entropy_from_prob(softmax_vec(eta))
    }

    summaries[[i]] <- data.frame(
      period_id = time_lookup$period_id[i],
      window_start = time_lookup$window_start[i],
      window_end = time_lookup$window_end[i],
      window_midpoint = time_lookup$window_midpoint[i],
      time_scaled = t_scaled,
      pe_mean = mean(pe),
      pe_q05 = unname(quantile(pe, 0.05)),
      pe_q50 = unname(quantile(pe, 0.50)),
      pe_q95 = unname(quantile(pe, 0.95))
    )
  }

  do.call(rbind, summaries)
}

summarise_window_pe <- function(fit, window_data) {
  pe_summary <- summarise_draw_matrix(as_fit_draw_matrix(fit, "normalized_pe"), "normalized_pe")
  eff_summary <- summarise_draw_matrix(
    as_fit_draw_matrix(fit, "effective_motif_n"),
    "effective_motif_n"
  )

  if (nrow(pe_summary) != nrow(window_data)) {
    stop(
      "The number of normalized_pe rows (", nrow(pe_summary),
      ") does not match window_data rows (", nrow(window_data), ")."
    )
  }

  out <- cbind(
    row_index = seq_len(nrow(window_data)),
    window_data,
    pe_model_mean = pe_summary$mean,
    pe_model_q05 = pe_summary$q05,
    pe_model_q50 = pe_summary$q50,
    pe_model_q95 = pe_summary$q95,
    effective_motif_n_mean = eff_summary$mean,
    effective_motif_n_q05 = eff_summary$q05,
    effective_motif_n_q50 = eff_summary$q50,
    effective_motif_n_q95 = eff_summary$q95
  )

  row.names(out) <- NULL
  out
}

summarise_pe_change <- function(fit, model = c("single", "multisite")) {
  model <- match.arg(model)
  vars <- if (model == "single") {
    c("normalized_pe_start", "normalized_pe_end", "normalized_pe_change")
  } else {
    c("population_pe_start", "population_pe_end", "population_pe_change")
  }

  draws <- as_fit_draw_matrix(fit, variables = vars)
  missing <- setdiff(vars, colnames(draws))
  if (length(missing)) {
    stop("Missing PE change columns: ", paste(missing, collapse = ", "))
  }

  data.frame(
    variable = vars,
    mean = colMeans(draws[, vars, drop = FALSE]),
    q05 = apply(draws[, vars, drop = FALSE], 2, quantile, probs = 0.05),
    q50 = apply(draws[, vars, drop = FALSE], 2, quantile, probs = 0.50),
    q95 = apply(draws[, vars, drop = FALSE], 2, quantile, probs = 0.95),
    row.names = NULL
  )
}

summarise_pe_fit <- function(fit, window_data, model = c("single", "multisite")) {
  model <- match.arg(model)

  list(
    window_pe = summarise_window_pe(fit, window_data),
    population_pe_time = summarise_population_pe_time(fit, window_data),
    pe_change = summarise_pe_change(fit, model = model)
  )
}
