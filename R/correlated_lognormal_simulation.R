simulate_correlated_lognormal_data <- function(
    n_pairs,
    sigma2,
    gamma,
    alpha,
    lambda,
    rho,
    beta0 = 1,
    beta1 = -0.01,
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

  frailty <- simulate_correlated_lognormal_frailty(
    pair_zygosity = pair_structure$pair_zygosity,
    sigma2 = sigma2,
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

  out$mean_w <- mean(
    frailty$w
  )

  out$var_w <- var(
    frailty$w
  )

  out$mean_frailty <- mean(
    frailty$frailty
  )

  out$var_frailty <- var(
    frailty$frailty
  )

  out$parameters <- list(
    sigma2 = sigma2,
    gamma = gamma,
    alpha = alpha,
    lambda = lambda,
    rho = rho,
    beta0 = beta0,
    beta1 = beta1
  )

  out
}
