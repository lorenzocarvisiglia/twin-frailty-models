simulate_power_gamma_data <- function(
    n_pairs,
    theta,
    gamma,
    alpha,
    lambda,
    rho,
    beta0 = 1,
    beta1 = -0.1,
    sd_pair = 0.1,
    sd_individual = 0.1,
    entry_max = 3,
    delayed_entry_prob = 0.5,
    censoring_max = 30,
    mz_proportion = 0.5,
    visit_step = 1
) {

  pair_structure <- make_pair_structure(
    n_pairs = n_pairs,
    mz_proportion = mz_proportion
  )

  frailty <- simulate_power_gamma_frailty(
    pair_zygosity = pair_structure$pair_zygosity,
    theta = theta,
    gamma = gamma
  )

  grid <- make_time_grid(
    followup_max = censoring_max,
    visit_step = visit_step
  )

  longitudinal <- simulate_longitudinal_covariate(
    n_pairs = n_pairs,
    grid_start = grid$start,
    beta0 = beta0,
    beta1 = beta1,
    sd_pair = sd_pair,
    sd_individual = sd_individual
  )

  event_time <- simulate_weibull_event_times(
    frailty = frailty$frailty,
    y = longitudinal$y,
    grid_start = grid$start,
    grid_stop = grid$stop,
    alpha = alpha,
    lambda = lambda,
    rho = rho
  )

  observation <- simulate_observation_process(
    event_time = event_time,
    pair_id = pair_structure$subjects$pair_original,
    delayed_entry_prob = delayed_entry_prob,
    entry_max = entry_max,
    censoring_max = censoring_max
  )

  out <- make_start_stop_data(
    pair_structure = pair_structure,
    y = longitudinal$y,
    grid_start = grid$start,
    grid_stop = grid$stop,
    event_time = event_time,
    observation = observation
  )

  keep <- observation$keep

  out$pair_zig <-
    pair_structure$pair_zygosity[
      observation$observed_pairs
    ]

  out$grid_start <- grid$start
  out$grid_stop <- grid$stop

  out$y_history <- longitudinal$y[
    keep,
    ,
    drop = FALSE
  ]

  out$n_subjects_observed <- sum(
    keep
  )

  out$n_events <- sum(
    observation$status[
      keep
    ]
  )

  out$event_proportion <-
    out$n_events /
    out$n_subjects_observed

  out$generated_frailty_mean <- mean(
    frailty$frailty
  )

  out$generated_frailty_variance <- var(
    frailty$frailty
  )

  out$parameters <- list(
    theta = theta,
    gamma = gamma,
    alpha = alpha,
    lambda = lambda,
    rho = rho,
    beta0 = beta0,
    beta1 = beta1
  )

  out
}


prepare_power_fit <- function(
    sim
) {

  prep <- prepare_frailty_fit(
    sim
  )

  required <- c(
    "grid_start",
    "grid_stop",
    "y_history"
  )

  missing <- setdiff(
    required,
    names(sim)
  )

  if (length(missing) > 0L) {
    stop(
      "simulation object is missing: ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }

  n_intervals <- length(
    sim$grid_start
  )

  if (length(sim$grid_stop) != n_intervals) {
    stop(
      "grid_start and grid_stop must have the same length"
    )
  }

  if (ncol(sim$y_history) != n_intervals) {
    stop(
      "y_history is incompatible with the time grid"
    )
  }

  if (nrow(sim$y_history) != prep$n_subjects) {
    stop(
      "y_history is incompatible with the observed subjects"
    )
  }

  entry_matrix <- matrix(
    prep$entry,
    nrow = prep$n_subjects,
    ncol = n_intervals
  )

  entry_interval_start <- matrix(
    sim$grid_start,
    nrow = prep$n_subjects,
    ncol = n_intervals,
    byrow = TRUE
  )

  entry_interval_stop <- matrix(
    sim$grid_stop,
    nrow = prep$n_subjects,
    ncol = n_intervals,
    byrow = TRUE
  )

  prep$y_history <- sim$y_history

  prep$entry_interval_start <-
    entry_interval_start

  prep$entry_interval_end <- pmin(
    entry_matrix,
    entry_interval_stop
  )

  prep
}
