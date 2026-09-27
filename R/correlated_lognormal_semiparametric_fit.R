make_correlated_semiparametric_quadrature <- function(
  n_nodes
) {
  rule <- normal_quadrature_rule(
    n_nodes
  )

  z <- rule$z
  log_weights <- rule$log_weights

  list(
    z = z,
    log_weights = log_weights,
    z1_2d = rep(
      z,
      each = n_nodes
    ),
    z2_2d = rep(
      z,
      times = n_nodes
    ),
    log_weights_2d =
      rep(
        log_weights,
        each = n_nodes
      ) +
      rep(
        log_weights,
        times = n_nodes
      ),
    n_nodes = n_nodes
  )
}


make_correlated_semiparametric_setup <- function(
  sim,
  n_nodes = 25L,
  reference_nodes = 60L
) {
  prep <-
    prepare_correlated_lognormal_fit(
      sim
    )

  event_rows <- which(
    prep$status == 1L
  )

  event_times <- sort(
    unique(
      prep$stop[
        event_rows
      ]
    )
  )

  if (length(event_times) == 0L) {
    stop("no observed event times")
  }

  event_index <- match(
    prep$stop[
      event_rows
    ],
    event_times
  )

  event_counts <- tabulate(
    event_index,
    nbins = length(
      event_times
    )
  )



  list(
    prep = prep,
    event_rows = event_rows,
    event_times = event_times,
    event_counts = event_counts,
    M = length(
      event_times
    ),
    rule =
      make_correlated_semiparametric_quadrature(
        n_nodes
      ),
    reference_rule =
      make_correlated_semiparametric_quadrature(
        reference_nodes
      ),
    n_nodes = n_nodes,
    reference_nodes =
      reference_nodes
  )
}


correlated_semiparametric_range_matrix <- function(
  subject_index,
  left,
  right,
  weight,
  n_subjects,
  M
) {
  z <- matrix(
    0,
    nrow = n_subjects,
    ncol = M + 1L
  )

  for (i in seq_along(subject_index)) {
    l <- left[i]
    r <- right[i]

    if (
      is.na(l) ||
      is.na(r) ||
      l > r ||
      r < 1L ||
      l > M
    ) {
      next
    }

    l <- max(
      1L,
      l
    )

    r <- min(
      M,
      r
    )

    s <- subject_index[i]
    w <- weight[i]

    z[s, l] <-
      z[s, l] +
      w

    if (r < M) {
      z[s, r + 1L] <-
        z[s, r + 1L] -
        w
    }
  }

  out <- z[
    ,
    seq_len(M),
    drop = FALSE
  ]

  for (s in seq_len(n_subjects)) {
    out[s, ] <- cumsum(
      out[s, ]
    )
  }

  out
}


correlated_semiparametric_weights <- function(
  setup,
  alpha
) {
  prep <- setup$prep
  event_times <- setup$event_times
  M <- setup$M

  entry_right <- findInterval(
    prep$entry,
    event_times
  )

  entry_weight <- exp(
    alpha *
      prep$y_entry
  )

  A_weight <-
    correlated_semiparametric_range_matrix(
      subject_index =
        seq_len(
          prep$n_subjects
        ),
      left =
        rep(
          1L,
          prep$n_subjects
        ),
      right =
        entry_right,
      weight =
        entry_weight,
      n_subjects =
        prep$n_subjects,
      M = M
    )

  post_left <-
    findInterval(
      prep$start,
      event_times
    ) +
    1L

  post_right <- findInterval(
    prep$stop,
    event_times
  )

  post_weight <- exp(
    alpha *
      prep$y
  )

  P_weight <-
    correlated_semiparametric_range_matrix(
      subject_index =
        prep$long_subject_index,
      left =
        post_left,
      right =
        post_right,
      weight =
        post_weight,
      n_subjects =
        prep$n_subjects,
      M = M
    )

  list(
    A = A_weight,
    R =
      A_weight +
      P_weight
  )
}


correlated_semiparametric_mz_integral <- function(
  d,
  H,
  sigma2,
  setup,
  reference = FALSE
) {
  d <- as.numeric(d)
  H <- as.numeric(H)

  if (
    length(d) != length(H) ||
    sigma2 < 0
  ) {
    return(
      rep(
        NA_real_,
        length(H)
      )
    )
  }

  if (sigma2 == 0) {
    return(
      -H
    )
  }

  rule <- if (reference) {
    setup$reference_rule
  } else {
    setup$rule
  }

  sigma <- sqrt(
    sigma2
  )

  w <- sigma *
    rule$z

  exp_w <- exp(
    w
  )

  log_integrand <-
    outer(
      d,
      w,
      "*"
    ) -
    outer(
      H,
      exp_w,
      "*"
    )

  log_integrand <- sweep(
    log_integrand,
    2,
    rule$log_weights,
    "+"
  )

  row_log_sum_exp(
    log_integrand
  )
}


correlated_semiparametric_dz_integral <- function(
  d1,
  d2,
  H1,
  H2,
  sigma2,
  gamma,
  setup,
  reference = FALSE
) {
  d1 <- as.numeric(d1)
  d2 <- as.numeric(d2)
  H1 <- as.numeric(H1)
  H2 <- as.numeric(H2)

  n <- length(d1)

  if (
    length(d2) != n ||
    length(H1) != n ||
    length(H2) != n ||
    sigma2 < 0 ||
    abs(gamma) > 1
  ) {
    return(
      rep(
        NA_real_,
        n
      )
    )
  }

  if (sigma2 == 0) {
    return(
      -H1 -
        H2
    )
  }

  if (gamma == 1) {
    return(
      correlated_semiparametric_mz_integral(
        d =
          d1 +
          d2,
        H =
          H1 +
          H2,
        sigma2 =
          sigma2,
        setup =
          setup,
        reference =
          reference
      )
    )
  }

  rule <- if (reference) {
    setup$reference_rule
  } else {
    setup$rule
  }

  sigma <- sqrt(
    sigma2
  )

  if (gamma == -1) {
    w1 <- sigma *
      rule$z

    w2 <- -w1

    exp_w1 <- exp(
      w1
    )

    exp_w2 <- exp(
      w2
    )

    log_integrand <-
      outer(
        d1,
        w1,
        "*"
      ) +
      outer(
        d2,
        w2,
        "*"
      ) -
      outer(
        H1,
        exp_w1,
        "*"
      ) -
      outer(
        H2,
        exp_w2,
        "*"
      )

    log_integrand <- sweep(
      log_integrand,
      2,
      rule$log_weights,
      "+"
    )

    return(
      row_log_sum_exp(
        log_integrand
      )
    )
  }

  w1 <- sigma *
    rule$z1_2d

  w2 <- sigma *
    (
      gamma *
        rule$z1_2d +
      sqrt(
        1 -
          gamma^2
      ) *
        rule$z2_2d
    )

  exp_w1 <- exp(
    w1
  )

  exp_w2 <- exp(
    w2
  )

  log_integrand <-
    outer(
      d1,
      w1,
      "*"
    ) +
    outer(
      d2,
      w2,
      "*"
    ) -
    outer(
      H1,
      exp_w1,
      "*"
    ) -
    outer(
      H2,
      exp_w2,
      "*"
    )

  log_integrand <- sweep(
    log_integrand,
    2,
    rule$log_weights_2d,
    "+"
  )

  row_log_sum_exp(
    log_integrand
  )
}


correlated_semiparametric_components <- function(
  setup,
  q,
  delta,
  reference = FALSE
) {
  prep <- setup$prep

  q <- as.numeric(q)
  delta <- as.numeric(delta)

  if (
    length(q) != 3L ||
    length(delta) != setup$M ||
    any(!is.finite(q)) ||
    any(!is.finite(delta)) ||
    any(delta <= 0)
  ) {
    return(
      list(
        loglik = -Inf
      )
    )
  }

  sigma2 <- q[1]
  gamma <- q[2]
  alpha <- q[3]

  if (
    sigma2 < 0 ||
    abs(gamma) > 1
  ) {
    return(
      list(
        loglik = -Inf
      )
    )
  }

  weights <-
    correlated_semiparametric_weights(
      setup = setup,
      alpha = alpha
    )

  A <- as.numeric(
    weights$A %*%
      delta
  )

  H <- as.numeric(
    weights$R %*%
      delta
  )

  if (
    any(!is.finite(A)) ||
    any(!is.finite(H)) ||
    any(A < 0) ||
    any(H < 0)
  ) {
    return(
      list(
        loglik = -Inf
      )
    )
  }

  log_event <-
    sum(
      setup$event_counts *
        log(delta)
    ) +
    alpha *
      sum(
        prep$y[
          setup$event_rows
        ]
      )

  loglik_mz <- 0

  if (
    length(
      prep$mz_index
    ) >
      0L
  ) {
    idx <- prep$mz_index

    s1 <- prep$subject1[
      idx
    ]

    s2 <- prep$subject2[
      idx
    ]

    d1 <- prep$subject_status[
      s1
    ]

    d2 <- prep$subject_status[
      s2
    ]

    d <- d1 +
      d2

    H_pair <-
      H[s1] +
      H[s2]

    A_pair <-
      A[s1] +
      A[s2]

    num <-
      correlated_semiparametric_mz_integral(
        d = d,
        H = H_pair,
        sigma2 = sigma2,
        setup = setup,
        reference = reference
      )

    den <-
      correlated_semiparametric_mz_integral(
        d =
          rep(
            0,
            length(idx)
          ),
        H = A_pair,
        sigma2 = sigma2,
        setup = setup,
        reference = reference
      )

    if (
      any(!is.finite(num)) ||
      any(!is.finite(den))
    ) {
      return(
        list(
          loglik = -Inf
        )
      )
    }

    loglik_mz <- sum(
      num -
        den
    )
  }

  loglik_dz <- 0

  if (
    length(
      prep$dz_index
    ) >
      0L
  ) {
    idx <- prep$dz_index

    s1 <- prep$subject1[
      idx
    ]

    s2 <- prep$subject2[
      idx
    ]

    d1 <- prep$subject_status[
      s1
    ]

    d2 <- prep$subject_status[
      s2
    ]

    num <-
      correlated_semiparametric_dz_integral(
        d1 = d1,
        d2 = d2,
        H1 = H[s1],
        H2 = H[s2],
        sigma2 = sigma2,
        gamma = gamma,
        setup = setup,
        reference = reference
      )

    den <-
      correlated_semiparametric_dz_integral(
        d1 =
          rep(
            0,
            length(idx)
          ),
        d2 =
          rep(
            0,
            length(idx)
          ),
        H1 = A[s1],
        H2 = A[s2],
        sigma2 = sigma2,
        gamma = gamma,
        setup = setup,
        reference = reference
      )

    if (
      any(!is.finite(num)) ||
      any(!is.finite(den))
    ) {
      return(
        list(
          loglik = -Inf
        )
      )
    }

    loglik_dz <- sum(
      num -
        den
    )
  }

  loglik <-
    log_event +
    loglik_mz +
    loglik_dz

  if (!is.finite(loglik)) {
    return(
      list(
        loglik = -Inf
      )
    )
  }

  list(
    loglik = loglik,
    log_event = log_event,
    loglik_mz = loglik_mz,
    loglik_dz = loglik_dz,
    A = A,
    H = H,
    weights = weights
  )
}


correlated_semiparametric_baseline_denominator <- function(
  setup,
  q,
  delta
) {
  prep <- setup$prep

  sigma2 <- q[1]
  gamma <- q[2]

  z <-
    correlated_semiparametric_components(
      setup = setup,
      q = q,
      delta = delta,
      reference = FALSE
    )

  if (!is.finite(z$loglik)) {
    return(
      rep(
        NA_real_,
        setup$M
      )
    )
  }

  A <- z$A
  H <- z$H

  A_weight <- z$weights$A
  R_weight <- z$weights$R

  denominator <- numeric(
    setup$M
  )

  if (
    length(
      prep$mz_index
    ) >
      0L
  ) {
    idx <- prep$mz_index

    s1 <- prep$subject1[
      idx
    ]

    s2 <- prep$subject2[
      idx
    ]

    d1 <- prep$subject_status[
      s1
    ]

    d2 <- prep$subject_status[
      s2
    ]

    d <- d1 +
      d2

    H_pair <-
      H[s1] +
      H[s2]

    A_pair <-
      A[s1] +
      A[s2]

    log_H <-
      correlated_semiparametric_mz_integral(
        d = d,
        H = H_pair,
        sigma2 = sigma2,
        setup = setup
      )

    log_H_next <-
      correlated_semiparametric_mz_integral(
        d = d + 1,
        H = H_pair,
        sigma2 = sigma2,
        setup = setup
      )

    zeros <- rep(
      0,
      length(idx)
    )

    ones <- rep(
      1,
      length(idx)
    )

    log_A <-
      correlated_semiparametric_mz_integral(
        d = zeros,
        H = A_pair,
        sigma2 = sigma2,
        setup = setup
      )

    log_A_next <-
      correlated_semiparametric_mz_integral(
        d = ones,
        H = A_pair,
        sigma2 = sigma2,
        setup = setup
      )

    ratio_H <- exp(
      log_H_next -
        log_H
    )

    ratio_A <- exp(
      log_A_next -
        log_A
    )

    if (
      any(!is.finite(ratio_H)) ||
      any(!is.finite(ratio_A))
    ) {
      return(
        rep(
          NA_real_,
          setup$M
        )
      )
    }

    R_pair <-
      R_weight[
        s1,
        ,
        drop = FALSE
      ] +
      R_weight[
        s2,
        ,
        drop = FALSE
      ]

    A_pair_weight <-
      A_weight[
        s1,
        ,
        drop = FALSE
      ] +
      A_weight[
        s2,
        ,
        drop = FALSE
      ]

    denominator <-
      denominator +
      as.numeric(
        crossprod(
          ratio_H,
          R_pair
        )
      ) -
      as.numeric(
        crossprod(
          ratio_A,
          A_pair_weight
        )
      )
  }

  if (
    length(
      prep$dz_index
    ) >
      0L
  ) {
    idx <- prep$dz_index

    s1 <- prep$subject1[
      idx
    ]

    s2 <- prep$subject2[
      idx
    ]

    d1 <- prep$subject_status[
      s1
    ]

    d2 <- prep$subject_status[
      s2
    ]

    log_H <-
      correlated_semiparametric_dz_integral(
        d1 = d1,
        d2 = d2,
        H1 = H[s1],
        H2 = H[s2],
        sigma2 = sigma2,
        gamma = gamma,
        setup = setup
      )

    log_H1 <-
      correlated_semiparametric_dz_integral(
        d1 = d1 + 1,
        d2 = d2,
        H1 = H[s1],
        H2 = H[s2],
        sigma2 = sigma2,
        gamma = gamma,
        setup = setup
      )

    log_H2 <-
      correlated_semiparametric_dz_integral(
        d1 = d1,
        d2 = d2 + 1,
        H1 = H[s1],
        H2 = H[s2],
        sigma2 = sigma2,
        gamma = gamma,
        setup = setup
      )

    zeros <- rep(
      0,
      length(idx)
    )

    ones <- rep(
      1,
      length(idx)
    )

    log_A <-
      correlated_semiparametric_dz_integral(
        d1 = zeros,
        d2 = zeros,
        H1 = A[s1],
        H2 = A[s2],
        sigma2 = sigma2,
        gamma = gamma,
        setup = setup
      )

    log_A1 <-
      correlated_semiparametric_dz_integral(
        d1 = ones,
        d2 = zeros,
        H1 = A[s1],
        H2 = A[s2],
        sigma2 = sigma2,
        gamma = gamma,
        setup = setup
      )

    log_A2 <-
      correlated_semiparametric_dz_integral(
        d1 = zeros,
        d2 = ones,
        H1 = A[s1],
        H2 = A[s2],
        sigma2 = sigma2,
        gamma = gamma,
        setup = setup
      )

    ratio_H1 <- exp(
      log_H1 -
        log_H
    )

    ratio_H2 <- exp(
      log_H2 -
        log_H
    )

    ratio_A1 <- exp(
      log_A1 -
        log_A
    )

    ratio_A2 <- exp(
      log_A2 -
        log_A
    )

    if (
      any(
        !is.finite(
          c(
            ratio_H1,
            ratio_H2,
            ratio_A1,
            ratio_A2
          )
        )
      )
    ) {
      return(
        rep(
          NA_real_,
          setup$M
        )
      )
    }

    denominator <-
      denominator +
      as.numeric(
        crossprod(
          ratio_H1,
          R_weight[
            s1,
            ,
            drop = FALSE
          ]
        )
      ) +
      as.numeric(
        crossprod(
          ratio_H2,
          R_weight[
            s2,
            ,
            drop = FALSE
          ]
        )
      ) -
      as.numeric(
        crossprod(
          ratio_A1,
          A_weight[
            s1,
            ,
            drop = FALSE
          ]
        )
      ) -
      as.numeric(
        crossprod(
          ratio_A2,
          A_weight[
            s2,
            ,
            drop = FALSE
          ]
        )
      )
  }

  denominator
}


fit_correlated_semiparametric_baseline <- function(
  setup,
  q,
  delta_start,
  max_iter = 200L,
  tolerance = 1e-6
) {
  delta <- as.numeric(
    delta_start
  )

  if (
    length(delta) != setup$M ||
    any(!is.finite(delta)) ||
    any(delta <= 0)
  ) {
    stop("invalid starting baseline")
  }

  current <-
    correlated_semiparametric_components(
      setup = setup,
      q = q,
      delta = delta
    )

  ll <- current$loglik

  if (!is.finite(ll)) {
    stop("invalid starting baseline likelihood")
  }

  converged <- FALSE
  max_score <- Inf
  iter <- 0L

  for (iter in seq_len(max_iter)) {
    denominator <-
      correlated_semiparametric_baseline_denominator(
        setup = setup,
        q = q,
        delta = delta
      )

    if (
      any(
        !is.finite(
          denominator
        )
      ) ||
      any(
        denominator <= 0
      )
    ) {
      stop("invalid baseline denominator")
    }

    target <-
      setup$event_counts /
      denominator

    if (
      any(
        !is.finite(
          target
        )
      ) ||
      any(
        target <= 0
      )
    ) {
      stop("invalid baseline target")
    }

    direction <-
      log(
        target
      ) -
      log(
        delta
      )

    omega <- 1
    accepted <- FALSE

    while (
      omega >=
        1 / 1024
    ) {
      candidate <- exp(
        log(
          delta
        ) +
        omega *
          direction
      )

      candidate_ll <-
        correlated_semiparametric_components(
          setup = setup,
          q = q,
          delta = candidate
        )$loglik

      if (
        is.finite(
          candidate_ll
        ) &&
        candidate_ll >=
          ll -
          1e-8
      ) {
        accepted <- TRUE
        break
      }

      omega <- omega /
        2
    }

    if (!accepted) {
      stop("baseline line search failed")
    }

    delta <- candidate
    ll <- candidate_ll

    denominator <-
      correlated_semiparametric_baseline_denominator(
        setup = setup,
        q = q,
        delta = delta
      )

    if (
      any(
        !is.finite(
          denominator
        )
      )
    ) {
      stop("invalid final baseline denominator")
    }

    score_eta <-
      setup$event_counts -
      delta *
        denominator

    max_score <- max(
      abs(
        score_eta
      )
    )

    if (
      max_score <
        tolerance
    ) {
      converged <- TRUE
      break
    }
  }

  list(
    delta = delta,
    loglik = ll,
    convergence =
      ifelse(
        converged,
        0L,
        1L
      ),
    iterations = iter,
    max_eta_score =
      max_score
  )
}

fit_correlated_lognormal_semiparametric_profile <- function(
  sim,
  verbose = FALSE
) {
  setup <-
    make_correlated_semiparametric_setup(
      sim
    )

  weibull_fit <-
    fit_correlated_lognormal(
      sim = sim
    )

  get_est <- function(
    name,
    fallback
  ) {
    z <- weibull_fit$estimates[[name]]

    if (
      length(z) == 1L &&
      is.finite(z)
    ) {
      as.numeric(z)
    } else {
      fallback
    }
  }

  sigma20 <- get_est(
    "sigma2",
    0.8
  )

  gamma0 <- get_est(
    "gamma",
    0.3
  )

  alpha0 <- get_est(
    "alpha",
    0
  )

  lambda0 <- get_est(
    "lambda",
    0.1
  )

  rho0 <- get_est(
    "rho",
    2
  )

  lower <- c(
    0,
    -1,
    -10
  )

  upper <- c(
    10,
    1,
    10
  )

  q0 <- c(
    sigma20,
    gamma0,
    alpha0
  )

  q0 <- pmin(
    upper,
    pmax(
      lower,
      q0
    )
  )

  previous_times <- c(
    0,
    setup$event_times[
      -setup$M
    ]
  )

  delta0 <-
    lambda0 *
      (
        setup$event_times^rho0 -
          previous_times^rho0
      )

  delta0 <- pmax(
    delta0,
    1e-10
  )

  structural_metric <- function(
    q,
    gradient
  ) {
    if (
      any(
        !is.finite(
          gradient
        )
      )
    ) {
      return(Inf)
    }

    sigma2_boundary <-
      q[1] -
        lower[1] <
      1e-6

    sigma2_upper_boundary <-
      upper[1] -
        q[1] <
      1e-6

    gamma_lower_boundary <-
      q[2] -
        lower[2] <
      1e-6

    gamma_upper_boundary <-
      upper[2] -
        q[2] <
      1e-6

    alpha_lower_boundary <-
      q[3] -
        lower[3] <
      1e-6

    alpha_upper_boundary <-
      upper[3] -
        q[3] <
      1e-6

    metric_components <- numeric()

    if (sigma2_boundary) {
      metric_components <- c(
        metric_components,
        max(
          -gradient[1],
          0
        )
      )
    } else if (sigma2_upper_boundary) {
      metric_components <- c(
        metric_components,
        max(
          gradient[1],
          0
        )
      )
    } else {
      metric_components <- c(
        metric_components,
        abs(
          gradient[1]
        )
      )
    }

    if (!sigma2_boundary) {
      if (gamma_lower_boundary) {
        metric_components <- c(
          metric_components,
          max(
            -gradient[2],
            0
          )
        )
      } else if (gamma_upper_boundary) {
        metric_components <- c(
          metric_components,
          max(
            gradient[2],
            0
          )
        )
      } else {
        metric_components <- c(
          metric_components,
          abs(
            gradient[2]
          )
        )
      }
    }

    if (alpha_lower_boundary) {
      metric_components <- c(
        metric_components,
        max(
          -gradient[3],
          0
        )
      )
    } else if (alpha_upper_boundary) {
      metric_components <- c(
        metric_components,
        max(
          gradient[3],
          0
        )
      )
    } else {
      metric_components <- c(
        metric_components,
        abs(
          gradient[3]
        )
      )
    }

    max(
      metric_components
    )
  }

  make_solver <- function(
    delta_initial,
    baseline_tolerance,
    baseline_max_iter,
    phase
  ) {
    state <- new.env(
      parent = emptyenv()
    )

    state$q <- NULL
    state$delta <-
      delta_initial
    state$profile <- NULL
    state$n_profile <- 0L
    state$history <- list()

    same_q <- function(
      x,
      y
    ) {
      !is.null(y) &&
        length(x) ==
          length(y) &&
        max(
          abs(
            x -
              y
          )
        ) <
          1e-13
    }

    profile_at <- function(
      q
    ) {
      q <- as.numeric(
        q
      )

      if (
        same_q(
          q,
          state$q
        )
      ) {
        return(
          state$profile
        )
      }

      baseline_fit <- tryCatch(
        fit_correlated_semiparametric_baseline(
          setup = setup,
          q = q,
          delta_start =
            state$delta,
          max_iter =
            baseline_max_iter,
          tolerance =
            baseline_tolerance
        ),
        error = function(e) {
          NULL
        }
      )

      if (
        is.null(
          baseline_fit
        ) ||
        !is.finite(
          baseline_fit$loglik
        )
      ) {
        baseline_fit <- tryCatch(
          fit_correlated_semiparametric_baseline(
            setup = setup,
            q = q,
            delta_start =
              delta_initial,
            max_iter =
              baseline_max_iter,
            tolerance =
              baseline_tolerance
          ),
          error = function(e) {
            NULL
          }
        )
      }

      if (
        is.null(
          baseline_fit
        ) ||
        !is.finite(
          baseline_fit$loglik
        )
      ) {
        return(NULL)
      }

      state$n_profile <-
        state$n_profile +
          1L

      out <- list(
        q = q,
        delta =
          baseline_fit$delta,
        loglik =
          baseline_fit$loglik,
        nll =
          -baseline_fit$loglik,
        baseline_iterations =
          baseline_fit$iterations,
        baseline_gradient =
          baseline_fit$max_eta_score,
        baseline_convergence =
          baseline_fit$convergence
      )

      state$q <- q
      state$delta <-
        baseline_fit$delta
      state$profile <- out

      state$history[[
        length(
          state$history
        ) +
          1L
      ]] <- list(
        phase = phase,
        evaluation =
          state$n_profile,
        sigma2 =
          q[1],
        gamma =
          q[2],
        alpha =
          q[3],
        logLik =
          baseline_fit$loglik,
        baseline_iterations =
          baseline_fit$iterations,
        baseline_gradient =
          baseline_fit$max_eta_score
      )

      if (verbose) {
        cat(
          phase,
          state$n_profile,
          "| sigma2",
          round(
            q[1],
            6
          ),
          "| gamma",
          round(
            q[2],
            6
          ),
          "| alpha",
          round(
            q[3],
            6
          ),
          "| logLik",
          round(
            baseline_fit$loglik,
            6
          ),
          "| baseline iter",
          baseline_fit$iterations,
          "| baseline score",
          signif(
            baseline_fit$max_eta_score,
            4
          ),
          "\n"
        )
      }

      out
    }

    fixed_delta_nll <- function(
      q,
      delta
    ) {
      z <-
        correlated_semiparametric_components(
          setup = setup,
          q = q,
          delta = delta
        )

      if (
        !is.finite(
          z$loglik
        )
      ) {
        return(
          1e100
        )
      }

      -z$loglik
    }

    envelope_gradient <- function(
      q,
      delta,
      rel_step = 1e-4
    ) {
      g <- numeric(
        length(q)
      )

      for (
        k in seq_along(q)
      ) {
        h <-
          rel_step *
            max(
              1,
              abs(
                q[k]
              )
            )

        q_minus <- q
        q_plus <- q

        q_minus[k] <-
          max(
            lower[k],
            q[k] -
              h
          )

        q_plus[k] <-
          min(
            upper[k],
            q[k] +
              h
          )

        if (
          q_plus[k] >
            q_minus[k]
        ) {
          f_minus <-
            fixed_delta_nll(
              q_minus,
              delta
            )

          f_plus <-
            fixed_delta_nll(
              q_plus,
              delta
            )

          if (
            is.finite(
              f_minus
            ) &&
            is.finite(
              f_plus
            )
          ) {
            g[k] <-
              (
                f_plus -
                  f_minus
              ) /
              (
                q_plus[k] -
                  q_minus[k]
              )
          } else {
            g[k] <- NA_real_
          }
        } else {
          g[k] <- 0
        }
      }

      g
    }

    list(
      profile_at =
        profile_at,
      envelope_gradient =
        envelope_gradient,
      state =
        state
    )
  }

  run_optimizer <- function(
    solver,
    q_start,
    factr,
    pgtol
  ) {
    initial <-
      solver$profile_at(
        q_start
      )

    if (
      is.null(
        initial
      )
    ) {
      stop(
        "initial profile failed"
      )
    }

    reference_nll <-
      initial$nll

    objective <- function(
      q
    ) {
      z <-
        solver$profile_at(
          q
        )

      if (
        is.null(
          z
        )
      ) {
        anchor <-
          if (
            is.null(
              solver$state$q
            )
          ) {
            q_start
          } else {
            solver$state$q
          }

        dq <- q -
          anchor

        return(
          1e6 +
            1e4 *
              sum(
                dq^2
              )
        )
      }

      z$nll -
        reference_nll
    }

    gradient <- function(
      q
    ) {
      z <-
        solver$profile_at(
          q
        )

      if (
        is.null(
          z
        )
      ) {
        anchor <-
          if (
            is.null(
              solver$state$q
            )
          ) {
            q_start
          } else {
            solver$state$q
          }

        return(
          2e4 *
            (
              q -
                anchor
            )
        )
      }

      solver$envelope_gradient(
        q = q,
        delta = z$delta,
        rel_step = 1e-4
      )
    }

    optim(
      par = q_start,
      fn = objective,
      gr = gradient,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(
        maxit = 200,
        factr = factr,
        pgtol = pgtol
      )
    )
  }

  main_solver <-
    make_solver(
      delta_initial =
        delta0,
      baseline_tolerance =
        1e-7,
      baseline_max_iter =
        1000L,
      phase =
        "profile"
    )

  main_fit <-
    run_optimizer(
      solver =
        main_solver,
      q_start =
        q0,
      factr =
        100,
      pgtol =
        1e-7
    )

  main_profile <-
    main_solver$profile_at(
      main_fit$par
    )

  if (
    is.null(
      main_profile
    )
  ) {
    stop(
      "final main profile failed"
    )
  }

  main_gradient <-
    main_solver$envelope_gradient(
      q =
        main_fit$par,
      delta =
        main_profile$delta,
      rel_step =
        1e-5
    )

  main_gradient_metric <-
    structural_metric(
      main_fit$par,
      main_gradient
    )

  polished <-
    main_gradient_metric >
      1e-4

  if (polished) {
    polish_solver <-
      make_solver(
        delta_initial =
          main_profile$delta,
        baseline_tolerance =
          1e-9,
        baseline_max_iter =
          2000L,
        phase =
          "polish"
      )

    polish_fit <-
      run_optimizer(
        solver =
          polish_solver,
        q_start =
          main_fit$par,
        factr =
          1,
        pgtol =
          1e-9
      )

    final_fit <-
      polish_fit

    final_solver <-
      polish_solver

    final_profile <-
      polish_solver$profile_at(
        polish_fit$par
      )

    history <-
      c(
        main_solver$state$history,
        polish_solver$state$history
      )

    profile_evaluations <-
      main_solver$state$n_profile +
      polish_solver$state$n_profile

    polishing_evaluations <-
      polish_solver$state$n_profile
  } else {
    final_fit <-
      main_fit

    final_solver <-
      main_solver

    final_profile <-
      main_profile

    history <-
      main_solver$state$history

    profile_evaluations <-
      main_solver$state$n_profile

    polishing_evaluations <-
      0L
  }

  if (
    is.null(
      final_profile
    )
  ) {
    stop(
      "final profile failed"
    )
  }

  q <-
    final_fit$par

  delta <-
    final_profile$delta

  final_gradient <-
    final_solver$envelope_gradient(
      q = q,
      delta = delta,
      rel_step = 1e-5
    )

  structural_gradient <-
    structural_metric(
      q,
      final_gradient
    )

  denominator <-
    correlated_semiparametric_baseline_denominator(
      setup = setup,
      q = q,
      delta = delta
    )

  final_baseline_score <-
    max(
      abs(
        setup$event_counts -
        delta *
          denominator
      )
    )

  final <-
    correlated_semiparametric_components(
      setup = setup,
      q = q,
      delta = delta
    )

  reference <-
    correlated_semiparametric_components(
      setup = setup,
      q = q,
      delta = delta,
      reference = TRUE
    )

  quadrature_difference <-
    final$loglik -
      reference$loglik

  sigma2_boundary <-
    q[1] <=
      1e-6

  sigma2_upper_boundary <-
    upper[1] -
      q[1] <
      1e-5

  gamma_lower_boundary <-
    q[2] -
      lower[2] <
      1e-5

  gamma_upper_boundary <-
    upper[2] -
      q[2] <
      1e-5

  gamma_boundary <-
    gamma_lower_boundary ||
    gamma_upper_boundary

  gamma_identifiable <-
    !sigma2_boundary

  alpha_lower_boundary <-
    q[3] -
      lower[3] <
      1e-5

  alpha_upper_boundary <-
    upper[3] -
      q[3] <
      1e-5

  boundary_problem <-
    sigma2_upper_boundary ||
    alpha_lower_boundary ||
    alpha_upper_boundary

  converged_outer <-
    main_fit$convergence ==
      0L

  converged_final <-
    final_fit$convergence ==
      0L &&
    is.finite(
      final$loglik
    ) &&
    structural_gradient <
      1e-4 &&
    final_baseline_score <
      1e-4

  fit_ok <-
    converged_final &&
    !boundary_problem &&
    is.finite(
      quadrature_difference
    ) &&
    abs(
      quadrature_difference
    ) <
      0.01

  list(
    fit_ok =
      fit_ok,
    converged_outer =
      converged_outer,
    converged_final =
      converged_final,
    final_iterations =
      profile_evaluations,
    sigma2 =
      q[1],
    gamma =
      q[2],
    alpha =
      q[3],
    delta =
      delta,
    event_times =
      setup$event_times,
    logLik =
      final$loglik,
    reference_logLik =
      reference$loglik,
    quadrature_difference =
      quadrature_difference,
    structural_gradient =
      structural_gradient,
    structural_gradient_vector =
      final_gradient,
    baseline_gradient =
      final_baseline_score,
    sigma2_boundary =
      sigma2_boundary,
    gamma_boundary =
      gamma_boundary,
    gamma_lower_boundary =
      gamma_lower_boundary,
    gamma_upper_boundary =
      gamma_upper_boundary,
    gamma_identifiable =
      gamma_identifiable,
    boundary_problem =
      boundary_problem,
    selected_method =
      if (
        polished
      ) {
        "profile_L_BFGS_B_polished"
      } else {
        "profile_L_BFGS_B"
      },
    history =
      history,
    weibull_initial_fit =
      weibull_fit,
    profile_evaluations =
      profile_evaluations,
    polishing_evaluations =
      polishing_evaluations,
    polished =
      polished,
    main_structural_gradient =
      main_gradient_metric,
    optimizer_convergence =
      final_fit$convergence,
    optimizer_message =
      final_fit$message
  )
}
