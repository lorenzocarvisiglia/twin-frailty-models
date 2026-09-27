make_pair_structure <- function(
    n_pairs,
    mz_proportion = 0.5
) {

  if (length(n_pairs) != 1 || n_pairs < 1) {
    stop("n_pairs must be a positive integer")
  }

  if (mz_proportion < 0 || mz_proportion > 1) {
    stop("mz_proportion must be between 0 and 1")
  }

  n_pairs <- as.integer(n_pairs)

  n_mz <- as.integer(
    round(n_pairs * mz_proportion)
  )

  n_mz <- max(
    0L,
    min(n_pairs, n_mz)
  )

  n_dz <- n_pairs - n_mz

  pair_zygosity <- c(
    rep("DZ", n_dz),
    rep("MZ", n_mz)
  )

  pair_id <- rep(
    seq_len(n_pairs),
    each = 2L
  )

  subject_id <- seq_len(
    2L * n_pairs
  )

  zig <- rep(
    pair_zygosity,
    each = 2L
  )

  subjects <- data.frame(
    pair_original = pair_id,
    id = subject_id,
    zig = zig,
    stringsAsFactors = FALSE
  )

  list(
    subjects = subjects,
    pair_zygosity = pair_zygosity,
    n_pairs = n_pairs,
    n_dz = n_dz,
    n_mz = n_mz
  )
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

  if (length(grid_start) < 1) {
    stop("grid_start must contain at least one time point")
  }

  if (sd_pair < 0 || sd_individual < 0) {
    stop("standard deviations must be non-negative")
  }

  n_pairs <- as.integer(n_pairs)
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


simulate_power_gamma_frailty <- function(
    pair_zygosity,
    theta,
    gamma
) {

  if (theta <= 0) {
    stop("theta must be positive")
  }

  if (gamma <= 0) {
    stop("gamma must be positive")
  }

  if (any(!pair_zygosity %in% c("DZ", "MZ"))) {
    stop("pair_zygosity must contain only DZ and MZ")
  }

  n_pairs <- length(pair_zygosity)

  v_pair <- rgamma(
    n_pairs,
    shape = 1 / theta,
    scale = theta
  )

  frailty_pair <- ifelse(
    pair_zygosity == "MZ",
    v_pair^gamma,
    v_pair
  )

  frailty <- rep(
    frailty_pair,
    each = 2L
  )

  list(
    frailty = frailty,
    frailty_pair = frailty_pair,
    v_pair = v_pair
  )
}


simulate_nested_gamma_frailty <- function(
    pair_zygosity,
    theta,
    phi
) {

  if (theta <= 0) {
    stop("theta must be positive")
  }

  if (phi <= 0) {
    stop("phi must be positive")
  }

  if (any(!pair_zygosity %in% c("DZ", "MZ"))) {
    stop("pair_zygosity must contain only DZ and MZ")
  }

  n_pairs <- length(pair_zygosity)
  n_subjects <- 2L * n_pairs

  v_pair <- rgamma(
    n_pairs,
    shape = 1 / theta,
    scale = theta
  )

  v <- rep(
    v_pair,
    each = 2L
  )

  w <- numeric(
    n_subjects
  )

  dz_pairs <- which(
    pair_zygosity == "DZ"
  )

  if (length(dz_pairs) > 0L) {

    dz_subjects <- as.vector(
      rbind(
        2L * dz_pairs - 1L,
        2L * dz_pairs
      )
    )

    w[dz_subjects] <- rgamma(
      length(dz_subjects),
      shape = 1 / phi,
      scale = phi
    )
  }

  mz_pairs <- which(
    pair_zygosity == "MZ"
  )

  if (length(mz_pairs) > 0L) {

    w_pair <- rgamma(
      length(mz_pairs),
      shape = 1 / phi,
      scale = phi
    )

    w[2L * mz_pairs - 1L] <- w_pair
    w[2L * mz_pairs] <- w_pair
  }

  frailty <- v * w

  list(
    frailty = frailty,
    v = v,
    v_pair = v_pair,
    w = w
  )
}


simulate_correlated_lognormal_frailty <- function(
    pair_zygosity,
    sigma2,
    gamma
) {

  if (sigma2 < 0) {
    stop("sigma2 must be non-negative")
  }

  if (abs(gamma) > 1) {
    stop("gamma must lie between -1 and 1")
  }

  if (any(!pair_zygosity %in% c("DZ", "MZ"))) {
    stop("pair_zygosity must contain only DZ and MZ")
  }

  n_pairs <- length(pair_zygosity)
  n_subjects <- 2L * n_pairs

  sigma <- sqrt(sigma2)

  w <- numeric(
    n_subjects
  )

  dz_pairs <- which(
    pair_zygosity == "DZ"
  )

  if (length(dz_pairs) > 0L) {

    z1 <- rnorm(
      length(dz_pairs)
    )

    z2 <- rnorm(
      length(dz_pairs)
    )

    w1 <- sigma * z1

    w2 <- sigma * (
      gamma * z1 +
        sqrt(1 - gamma^2) * z2
    )

    w[2L * dz_pairs - 1L] <- w1
    w[2L * dz_pairs] <- w2
  }

  mz_pairs <- which(
    pair_zygosity == "MZ"
  )

  if (length(mz_pairs) > 0L) {

    w_pair <- rnorm(
      length(mz_pairs),
      mean = 0,
      sd = sigma
    )

    w[2L * mz_pairs - 1L] <- w_pair
    w[2L * mz_pairs] <- w_pair
  }

  frailty <- exp(w)

  list(
    frailty = frailty,
    w = w
  )
}


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

  if (any(frailty <= 0)) {
    stop("frailty values must be positive")
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

  event_time <- rep(
    Inf,
    n
  )

  active <- rep(
    TRUE,
    n
  )

  for (k in seq_along(grid_start)) {

    if (!any(active)) {
      break
    }

    a <- grid_start[k]
    b <- grid_stop[k]

    multiplier <- frailty *
      lambda *
      exp(
        alpha * y[, k]
      )

    increment <- multiplier *
      (
        b^rho -
          a^rho
      )

    hit <- active &
      target <=
      cumulative +
      increment

    if (any(hit)) {

      remaining <- target[hit] -
        cumulative[hit]

      event_time[hit] <- (
        a^rho +
          remaining /
          multiplier[hit]
      )^(1 / rho)

      active[hit] <- FALSE
    }

    cumulative[active] <-
      cumulative[active] +
      increment[active]
  }

  event_time
}


simulate_observation_process <- function(
    event_time,
    pair_id,
    delayed_entry_prob = 0.5,
    entry_max = 3,
    censoring_max = 30
) {

  if (length(event_time) != length(pair_id)) {
    stop("event_time and pair_id must have the same length")
  }

  if (
    delayed_entry_prob < 0 ||
    delayed_entry_prob > 1
  ) {
    stop("delayed_entry_prob must be between 0 and 1")
  }

  if (entry_max < 0) {
    stop("entry_max must be non-negative")
  }

  if (censoring_max <= 0) {
    stop("censoring_max must be positive")
  }

  pair_sizes <- table(pair_id)

  if (any(pair_sizes != 2L)) {
    stop("each pair must contain exactly two individuals")
  }

  n_subjects <- length(event_time)

  delayed <- rbinom(
    n_subjects,
    size = 1,
    prob = delayed_entry_prob
  )

  entry <- ifelse(
    delayed == 1L,
    runif(
      n_subjects,
      min = 0,
      max = entry_max
    ),
    0
  )

  censoring <- runif(
    n_subjects,
    min = 0,
    max = censoring_max
  )

  observed_stop <- pmin(
    event_time,
    censoring
  )

  status <- as.integer(
    event_time <= censoring
  )

  eligible <- entry < observed_stop

  pair_levels <- unique(pair_id)

  pair_is_observed <- vapply(
    pair_levels,
    function(p) {
      all(
        eligible[pair_id == p]
      )
    },
    logical(1)
  )

  observed_pairs <- pair_levels[
    pair_is_observed
  ]

  keep <- pair_id %in%
    observed_pairs

  list(
    entry = entry,
    censoring = censoring,
    observed_stop = observed_stop,
    status = status,
    eligible = eligible,
    keep = keep,
    observed_pairs = observed_pairs
  )
}


make_start_stop_data <- function(
    pair_structure,
    y,
    grid_start,
    grid_stop,
    event_time,
    observation
) {

  subjects_all <- pair_structure$subjects

  if (nrow(subjects_all) != nrow(y)) {
    stop("pair_structure and y contain different numbers of individuals")
  }

  if (length(grid_start) != ncol(y)) {
    stop("the time grid is incompatible with y")
  }

  if (length(grid_stop) != ncol(y)) {
    stop("the time grid is incompatible with y")
  }

  keep_index <- which(
    observation$keep
  )

  if (length(keep_index) == 0L) {
    stop("no complete pairs remain after delayed-entry filtering")
  }

  subjects_keep <- subjects_all[
    keep_index,
    ,
    drop = FALSE
  ]

  observed_pairs <- observation$observed_pairs

  pair_index <- match(
    subjects_keep$pair_original,
    observed_pairs
  )

  entry_keep <- observation$entry[
    keep_index
  ]

  stop_keep <- observation$observed_stop[
    keep_index
  ]

  status_keep <- observation$status[
    keep_index
  ]

  event_time_keep <- event_time[
    keep_index
  ]

  censoring_keep <- observation$censoring[
    keep_index
  ]

  y_keep <- y[
    keep_index,
    ,
    drop = FALSE
  ]

  n_keep <- length(
    keep_index
  )

  n_intervals <- length(
    grid_start
  )

  entry_matrix <- matrix(
    entry_keep,
    nrow = n_keep,
    ncol = n_intervals
  )

  stop_matrix <- matrix(
    stop_keep,
    nrow = n_keep,
    ncol = n_intervals
  )

  grid_start_matrix <- matrix(
    grid_start,
    nrow = n_keep,
    ncol = n_intervals,
    byrow = TRUE
  )

  grid_stop_matrix <- matrix(
    grid_stop,
    nrow = n_keep,
    ncol = n_intervals,
    byrow = TRUE
  )

  interval_start <- pmax(
    entry_matrix,
    grid_start_matrix
  )

  interval_stop <- pmin(
    stop_matrix,
    grid_stop_matrix
  )

  mask <- interval_start <
    interval_stop

  rr <- row(
    interval_start
  )[mask]

  cc <- col(
    interval_start
  )[mask]

  start_long <- interval_start[
    mask
  ]

  stop_long <- interval_stop[
    mask
  ]

  status_long <- as.integer(
    status_keep[rr] == 1L &
      abs(
        stop_long -
          stop_keep[rr]
      ) <
      sqrt(.Machine$double.eps)
  )

  long <- data.frame(
    pair_index = pair_index[rr],
    pair_original = subjects_keep$pair_original[rr],
    id = subjects_keep$id[rr],
    zig = subjects_keep$zig[rr],
    start = start_long,
    stop = stop_long,
    status = status_long,
    y = y_keep[
      cbind(
        rr,
        cc
      )
    ],
    stringsAsFactors = FALSE
  )

  subjects <- data.frame(
    pair_index = pair_index,
    pair_original = subjects_keep$pair_original,
    id = subjects_keep$id,
    zig = subjects_keep$zig,
    entry = entry_keep,
    stop = stop_keep,
    status = status_keep,
    event_time = event_time_keep,
    censoring = censoring_keep,
    stringsAsFactors = FALSE
  )

  list(
    long = long,
    subjects = subjects,
    n_pairs_generated = pair_structure$n_pairs,
    n_pairs_observed = length(observed_pairs),
    n_dz_observed = sum(
      pair_structure$pair_zygosity[
        observed_pairs
      ] == "DZ"
    ),
    n_mz_observed = sum(
      pair_structure$pair_zygosity[
        observed_pairs
      ] == "MZ"
    )
  )
}
