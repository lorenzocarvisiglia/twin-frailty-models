simulate_weibull_event_times <- function(
    frailty,
    y,
    grid_start,
    grid_stop,
    alpha,
    lambda,
    rho
) {

  if (lambda <= 0 || rho <= 0) {
    stop("lambda and rho must be positive")
  }

  if (length(grid_start) != length(grid_stop)) {
    stop("grid_start and grid_stop must have the same length")
  }

  if (nrow(y) != length(frailty)) {
    stop("the number of rows of y must equal the length of frailty")
  }

  if (ncol(y) != length(grid_start)) {
    stop("the number of columns of y must equal the number of time intervals")
  }

  n <- length(frailty)

  target <- rexp(n)
  cumulative <- numeric(n)
  event_time <- rep(Inf, n)
  active <- rep(TRUE, n)

  for (k in seq_along(grid_start)) {

    if (!any(active)) {
      break
    }

    a <- grid_start[k]
    b <- grid_stop[k]

    multiplier <- frailty * lambda * exp(alpha * y[, k])

    increment <- multiplier * (b^rho - a^rho)

    hit <- active & target <= cumulative + increment

    if (any(hit)) {

      remaining <- target[hit] - cumulative[hit]

      event_time[hit] <- (
        a^rho +
          remaining / multiplier[hit]
      )^(1 / rho)

      active[hit] <- FALSE
    }

    cumulative[active] <-
      cumulative[active] +
      increment[active]
  }

  event_time
}


make_time_grid <- function(
    followup_max,
    visit_step = 1
) {

  if (followup_max <= 0) {
    stop("followup_max must be positive")
  }

  if (visit_step <= 0) {
    stop("visit_step must be positive")
  }

  grid_start <- seq(
    0,
    followup_max - 1e-12,
    by = visit_step
  )

  grid_stop <- pmin(
    grid_start + visit_step,
    followup_max
  )

  list(
    start = grid_start,
    stop = grid_stop
  )
}


simulate_longitudinal_covariate <- function(
    n_pairs,
    grid_start,
    beta0,
    beta1,
    sd_pair = 0.1,
    sd_individual = 0.1
) {

  if (n_pairs <= 0) {
    stop("n_pairs must be positive")
  }

  if (sd_pair < 0 || sd_individual < 0) {
    stop("standard deviations must be non-negative")
  }

  n_subjects <- 2L * n_pairs

  pair_effect <- rnorm(
    n_pairs,
    mean = 0,
    sd = sd_pair
  )

  individual_effect <- rnorm(
    n_subjects,
    mean = 0,
    sd = sd_individual
  )

  intercept <- beta0 +
    rep(pair_effect, each = 2L) +
    individual_effect

  y <- matrix(
    intercept,
    nrow = n_subjects,
    ncol = length(grid_start)
  ) +
    matrix(
      beta1 * grid_start,
      nrow = n_subjects,
      ncol = length(grid_start),
      byrow = TRUE
    )

  list(
    y = y,
    pair_effect = pair_effect,
    individual_effect = individual_effect
  )
}

