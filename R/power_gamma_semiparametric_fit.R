make_power_semiparametric_setup <- function(
  sim,
  n_nodes = 64L
) {
  prep <- prepare_power_fit(
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
    M = length(event_times),
    rule = power_legendre_rule(
      n_nodes
    ),
    n_nodes = n_nodes
  )
}

power_semiparametric_range_matrix <- function(
  pair_index,
  left,
  right,
  weight,
  n_pairs,
  M
) {
  z <- matrix(
    0,
    nrow = n_pairs,
    ncol = M + 1L
  )

  for (i in seq_along(pair_index)) {
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

    p <- pair_index[i]
    w <- weight[i]

    z[p, l] <-
      z[p, l] +
      w

    if (r < M) {
      z[p, r + 1L] <-
        z[p, r + 1L] -
        w
    }
  }

  out <- z[
    ,
    seq_len(M),
    drop = FALSE
  ]

  for (p in seq_len(n_pairs)) {
    out[p, ] <- cumsum(
      out[p, ]
    )
  }

  out
}

power_semiparametric_weights <- function(
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
    power_semiparametric_range_matrix(
      pair_index =
        prep$pair_index_subject,
      left =
        rep(
          1L,
          length(
            prep$entry
          )
        ),
      right =
        entry_right,
      weight =
        entry_weight,
      n_pairs =
        prep$n_pairs,
      M =
        M
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
    power_semiparametric_range_matrix(
      pair_index =
        prep$pair_index_long,
      left =
        post_left,
      right =
        post_right,
      weight =
        post_weight,
      n_pairs =
        prep$n_pairs,
      M =
        M
    )

  list(
    A = A_weight,
    R = A_weight + P_weight
  )
}

power_semiparametric_components <- function(
  setup,
  q,
  delta,
  reference = FALSE
) {
  prep <- setup$prep

  theta <- exp(
    q[1]
  )

  gamma <- exp(
    q[2]
  )

  alpha <- q[3]

  weights <-
    power_semiparametric_weights(
      setup,
      alpha
    )

  A <- as.numeric(
    weights$A %*%
      delta
  )

  H <- as.numeric(
    weights$R %*%
      delta
  )

  log_event <-
    sum(
      setup$event_counts *
        log(
          delta
        )
    ) +
    alpha *
      sum(
        prep$y[
          setup$event_rows
        ]
      )

  shape <- 1 /
    theta

  loglik_dz <- 0

  if (
    length(
      prep$dz_index
    ) >
      0L
  ) {
    idx <- prep$dz_index

    d <- prep$d_pair[
      idx
    ]

    loglik_dz <- sum(
      d *
        log(theta) +
        lgamma(
          shape +
            d
        ) -
        lgamma(
          shape
        ) +
        shape *
          log1p(
            theta *
              A[idx]
          ) -
        (
          shape +
            d
        ) *
          log1p(
            theta *
              H[idx]
          )
    )
  }

  loglik_mz <- 0

  if (
    length(
      prep$mz_index
    ) >
      0L
  ) {
    idx <- prep$mz_index

    d <- prep$d_pair[
      idx
    ]

    if (!reference) {
      num <-
        power_mz_log_integral_asymmetric(
          B = H[idx],
          d = d,
          theta = theta,
          gamma = gamma,
          rule = setup$rule,
          log_drop = 20
        )

      den <-
        power_mz_log_integral_asymmetric(
          B = A[idx],
          d = rep(
            0,
            length(idx)
          ),
          theta = theta,
          gamma = gamma,
          rule = setup$rule,
          log_drop = 20
        )
    } else {
      num <- vapply(
        seq_along(idx),
        function(k) {
          power_mz_log_integral_reference(
            B = H[idx[k]],
            d = d[k],
            theta = theta,
            gamma = gamma
          )
        },
        numeric(1)
      )

      den <- vapply(
        seq_along(idx),
        function(k) {
          power_mz_log_integral_reference(
            B = A[idx[k]],
            d = 0,
            theta = theta,
            gamma = gamma
          )
        },
        numeric(1)
      )
    }

    if (
      any(
        !is.finite(
          num
        )
      ) ||
      any(
        !is.finite(
          den
        )
      )
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

  list(
    loglik =
      log_event +
        loglik_dz +
        loglik_mz,
    A = A,
    H = H,
    weights = weights,
    loglik_event =
      log_event,
    loglik_dz =
      loglik_dz,
    loglik_mz =
      loglik_mz
  )
}

power_semiparametric_baseline_denominator <- function(
  setup,
  q,
  delta
) {
  prep <- setup$prep

  theta <- exp(
    q[1]
  )

  gamma <- exp(
    q[2]
  )

  alpha <- q[3]

  z <- power_semiparametric_components(
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
      prep$dz_index
    ) >
      0L
  ) {
    idx <- prep$dz_index

    d <- prep$d_pair[
      idx
    ]

    cH <-
      (
        1 +
          theta *
            d
      ) /
      (
        1 +
          theta *
            H[idx]
      )

    cA <-
      1 /
      (
        1 +
          theta *
            A[idx]
      )

    denominator <-
      denominator +
      as.numeric(
        crossprod(
          cH,
          R_weight[
            idx,
            ,
            drop = FALSE
          ]
        )
      ) -
      as.numeric(
        crossprod(
          cA,
          A_weight[
            idx,
            ,
            drop = FALSE
          ]
        )
      )
  }

  if (
    length(
      prep$mz_index
    ) >
      0L
  ) {
    idx <- prep$mz_index

    d <- prep$d_pair[
      idx
    ]

    log_I_H_d <-
      power_mz_log_integral_asymmetric(
        B = H[idx],
        d = d,
        theta = theta,
        gamma = gamma,
        rule = setup$rule,
        log_drop = 20
      )

    log_I_H_next <-
      power_mz_log_integral_asymmetric(
        B = H[idx],
        d = d + 1L,
        theta = theta,
        gamma = gamma,
        rule = setup$rule,
        log_drop = 20
      )

    log_I_A_0 <-
      power_mz_log_integral_asymmetric(
        B = A[idx],
        d = rep(
          0,
          length(idx)
        ),
        theta = theta,
        gamma = gamma,
        rule = setup$rule,
        log_drop = 20
      )

    log_I_A_1 <-
      power_mz_log_integral_asymmetric(
        B = A[idx],
        d = rep(
          1,
          length(idx)
        ),
        theta = theta,
        gamma = gamma,
        rule = setup$rule,
        log_drop = 20
      )

    ratio_H <- exp(
      log_I_H_next -
        log_I_H_d
    )

    ratio_A <- exp(
      log_I_A_1 -
        log_I_A_0
    )

    denominator <-
      denominator +
      as.numeric(
        crossprod(
          ratio_H,
          R_weight[
            idx,
            ,
            drop = FALSE
          ]
        )
      ) -
      as.numeric(
        crossprod(
          ratio_A,
          A_weight[
            idx,
            ,
            drop = FALSE
          ]
        )
      )
  }

  denominator
}

fit_power_semiparametric_baseline <- function(
  setup,
  q,
  delta_start,
  max_iter = 200L,
  tolerance = 1e-6
) {
  delta <- delta_start

  current <-
    power_semiparametric_components(
      setup = setup,
      q = q,
      delta = delta
    )

  ll <- current$loglik

  if (!is.finite(ll)) {
    stop(
      "invalid starting baseline likelihood"
    )
  }

  converged <- FALSE
  max_score <- Inf

  for (iter in seq_len(max_iter)) {
    denominator <-
      power_semiparametric_baseline_denominator(
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
      stop(
        "invalid baseline denominator"
      )
    }

    target <-
      setup$event_counts /
        denominator

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
        power_semiparametric_components(
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
      stop(
        "baseline line search failed"
      )
    }

    delta <- candidate
    ll <- candidate_ll

    denominator <-
      power_semiparametric_baseline_denominator(
        setup = setup,
        q = q,
        delta = delta
      )

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

optimize_power_semiparametric_structural <- function(
  setup,
  delta,
  q_start
) {
  lower <- c(
    log(0.02),
    log(0.05),
    -10
  )

  upper <- c(
    log(10),
    log(8),
    10
  )

  nll <- function(q) {
    z <-
      power_semiparametric_components(
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

  fit1 <- tryCatch(
    nlminb(
      start = q_start,
      objective = nll,
      lower = lower,
      upper = upper,
      control = list(
        eval.max = 1500,
        iter.max = 1000,
        rel.tol = 1e-9,
        x.tol = 1e-9
      )
    ),
    error = function(e) {
      NULL
    }
  )

  candidates <- list()

  if (
    !is.null(
      fit1
    ) &&
    is.finite(
      fit1$objective
    )
  ) {
    candidates$nlminb <- list(
      par = fit1$par,
      objective = fit1$objective,
      convergence = fit1$convergence,
      message = fit1$message
    )

    nll_reference <-
      nll(
        fit1$par
      )

    nll_centered <- function(q) {
      nll(q) -
        nll_reference
    }

    fit2 <- tryCatch(
      optim(
        par = fit1$par,
        fn = nll_centered,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(
          maxit = 5000,
          factr = 10,
          pgtol = 1e-10,
          ndeps = rep(
            1e-4,
            3L
          )
        )
      ),
      error = function(e) {
        NULL
      }
    )

    if (
      !is.null(
        fit2
      ) &&
      all(
        is.finite(
          fit2$par
        )
      )
    ) {
      fit2_objective <-
        nll(
          fit2$par
        )

      if (
        is.finite(
          fit2_objective
        )
      ) {
        candidates$L_BFGS_B <- list(
          par = fit2$par,
          objective = fit2_objective,
          convergence = fit2$convergence,
          message = fit2$message
        )
      }
    }
  }

  if (length(candidates) == 0L) {
    stop(
      "structural optimization failed"
    )
  }

  ok <- vapply(
    candidates,
    function(x) {
      x$convergence == 0L
    },
    logical(1)
  )

  objective <- vapply(
    candidates,
    function(x) {
      x$objective
    },
    numeric(1)
  )

  candidate_gradient <- lapply(
    candidates,
    function(x) {
      raw <- tryCatch(
        bounded_gradient(
          fn = nll,
          x = x$par,
          lower = lower,
          upper = upper,
          rel_step = 1e-4
        ),
        error = function(e) {
          rep(
            NA_real_,
            3L
          )
        }
      )

      at_lower <-
        x$par <=
          lower + 1e-5

      at_upper <-
        x$par >=
          upper - 1e-5

      projected <- raw

      projected[
        at_lower &
          raw > 0
      ] <- 0

      projected[
        at_upper &
          raw < 0
      ] <- 0

      projected
    }
  )

  candidate_max_gradient <- vapply(
    candidate_gradient,
    function(g) {
      if (
        all(
          is.finite(
            g
          )
        )
      ) {
        max(
          abs(
            g
          )
        )
      } else {
        Inf
      }
    },
    numeric(1)
  )

  pool <- if (any(ok)) {
    which(ok)
  } else {
    seq_along(
      candidates
    )
  }

  best_objective <- min(
    objective[
      pool
    ]
  )

  objective_tolerance <- 1e-6

  near_best <- pool[
    objective[
      pool
    ] <=
      best_objective +
        objective_tolerance
  ]

  if (
    length(
      near_best
    ) > 1L &&
      any(
        is.finite(
          candidate_max_gradient[
            near_best
          ]
        )
      )
  ) {
    selected <- near_best[
      which.min(
        candidate_max_gradient[
          near_best
        ]
      )
    ]
  } else {
    selected <- pool[
      which.min(
        objective[
          pool
        ]
      )
    ]
  }

  fit <- candidates[[
    selected
  ]]

  q_hat <- fit$par

  raw_gradient <- tryCatch(
    bounded_gradient(
      fn = nll,
      x = q_hat,
      lower = lower,
      upper = upper,
      rel_step = 1e-5
    ),
    error = function(e) {
      rep(
        NA_real_,
        3L
      )
    }
  )

  at_lower <-
    q_hat <=
      lower + 1e-5

  at_upper <-
    q_hat >=
      upper - 1e-5

  gradient <- raw_gradient

  gradient[
    at_lower &
      raw_gradient > 0
  ] <- 0

  gradient[
    at_upper &
      raw_gradient < 0
  ] <- 0

  max_abs_gradient <- if (
    all(
      is.finite(
        gradient
      )
    )
  ) {
    max(
      abs(
        gradient
      )
    )
  } else {
    Inf
  }

  list(
    par = q_hat,
    objective = fit$objective,
    convergence = fit$convergence,
    message = fit$message,
    selected_method =
      names(candidates)[selected],
    max_abs_gradient =
      max_abs_gradient,
    gradient =
      gradient
  )
}

fit_power_gamma_semiparametric <- function(
  sim,
  outer_max_iter = 100L,
  outer_tolerance = 1e-5,
  verbose = FALSE
) {
  setup <-
    make_power_semiparametric_setup(
      sim
    )

  weibull_fit <- fit_power_gamma(
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

  theta0 <- get_est(
    "theta",
    0.5
  )

  gamma0 <- get_est(
    "gamma",
    1
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

  q <- c(
    log(theta0),
    log(gamma0),
    alpha0
  )

  previous_times <- c(
    0,
    setup$event_times[
      -setup$M
    ]
  )

  delta <-
    lambda0 *
      (
        setup$event_times^rho0 -
          previous_times^rho0
      )

  delta <- pmax(
    delta,
    1e-10
  )

  loglik_previous <- -Inf
  converged <- FALSE

  history <- vector(
    "list",
    outer_max_iter
  )

  for (
    outer_iter in
      seq_len(
        outer_max_iter
      )
  ) {
    q_old <- q
    delta_old <- delta

    baseline_fit <-
      fit_power_semiparametric_baseline(
        setup = setup,
        q = q,
        delta_start = delta,
        max_iter = 200L,
        tolerance = 1e-6
      )

    delta <- baseline_fit$delta

    structural_fit <-
      optimize_power_semiparametric_structural(
        setup = setup,
        delta = delta,
        q_start = q
      )

    q <- structural_fit$par

    current <-
      power_semiparametric_components(
        setup = setup,
        q = q,
        delta = delta
      )

    denominator <-
      power_semiparametric_baseline_denominator(
        setup = setup,
        q = q,
        delta = delta
      )

    baseline_score <-
      setup$event_counts -
        delta *
          denominator

    max_baseline_score <- max(
      abs(
        baseline_score
      )
    )

    q_change <- max(
      abs(
        q -
          q_old
      )
    )

    delta_change <- max(
      abs(
        log(
          delta
        ) -
          log(
            delta_old
          )
      )
    )

    loglik_change <- if (
      is.finite(
        loglik_previous
      )
    ) {
      abs(
        current$loglik -
          loglik_previous
      )
    } else {
      Inf
    }

    history[[outer_iter]] <-
      data.frame(
        iteration =
          outer_iter,
        logLik =
          current$loglik,
        theta =
          exp(
            q[1]
          ),
        gamma =
          exp(
            q[2]
          ),
        alpha =
          q[3],
        q_change =
          q_change,
        delta_change =
          delta_change,
        structural_gradient =
          structural_fit$max_abs_gradient,
        baseline_gradient =
          max_baseline_score,
        baseline_iterations =
          baseline_fit$iterations
      )

    if (verbose) {
      cat(
        "outer",
        outer_iter,
        "| logLik",
        format(
          current$loglik,
          digits = 12
        ),
        "| theta",
        format(
          exp(
            q[1]
          ),
          digits = 7
        ),
        "| gamma",
        format(
          exp(
            q[2]
          ),
          digits = 7
        ),
        "| alpha",
        format(
          q[3],
          digits = 7
        ),
        "| q change",
        format(
          q_change,
          scientific = TRUE,
          digits = 4
        ),
        "| baseline score",
        format(
          max_baseline_score,
          scientific = TRUE,
          digits = 4
        ),
        "| structural score",
        format(
          structural_fit$max_abs_gradient,
          scientific = TRUE,
          digits = 4
        ),
        "\n"
      )
    }

    if (
      outer_iter > 1L &&
      q_change <
        outer_tolerance &&
      max_baseline_score <
        1e-4 &&
      structural_fit$max_abs_gradient <
        1e-4 &&
      loglik_change <
        1e-5
    ) {
      converged <- TRUE
      break
    }

    loglik_previous <-
      current$loglik
  }

  history <- do.call(
    rbind,
    history[
      seq_len(
        outer_iter
      )
    ]
  )

  final_converged <- FALSE
  final_max_iter <- 50L
  final_loglik_previous <- -Inf

  for (
    final_iter in
      seq_len(
        final_max_iter
      )
  ) {
    q_old_final <- q

    final_baseline <-
      fit_power_semiparametric_baseline(
        setup = setup,
        q = q,
        delta_start = delta,
        max_iter = 300L,
        tolerance = 1e-7
      )

    delta <-
      final_baseline$delta

    final_structural <-
      optimize_power_semiparametric_structural(
        setup = setup,
        delta = delta,
        q_start = q
      )

    q <- final_structural$par

    final <-
      power_semiparametric_components(
        setup = setup,
        q = q,
        delta = delta
      )

    denominator <-
      power_semiparametric_baseline_denominator(
        setup = setup,
        q = q,
        delta = delta
      )

    final_baseline_score <- max(
      abs(
        setup$event_counts -
          delta *
            denominator
      )
    )

    final_q_change <- max(
      abs(
        q -
          q_old_final
      )
    )

    final_loglik_change <- if (
      is.finite(
        final_loglik_previous
      )
    ) {
      abs(
        final$loglik -
          final_loglik_previous
      )
    } else {
      Inf
    }

    if (
      final_iter > 1L &&
      final_q_change <
        outer_tolerance &&
      final_baseline_score <
        1e-4 &&
      final_structural$max_abs_gradient <
        1e-4 &&
      final_loglik_change <
        1e-5
    ) {
      final_converged <- TRUE
      break
    }

    final_loglik_previous <-
      final$loglik
  }

  reference <-
    power_semiparametric_components(
      setup = setup,
      q = q,
      delta = delta,
      reference = TRUE
    )

  quadrature_difference <-
    final$loglik -
      reference$loglik

  lower <- c(
    log(0.02),
    log(0.05),
    -10
  )

  upper <- c(
    log(10),
    log(8),
    10
  )

  boundary_hit <- any(
    q -
      lower <
      1e-5 |
      upper -
        q <
        1e-5
  )

  fit_ok <-
    final_converged &&
    final_structural$convergence ==
      0L &&
    !boundary_hit &&
    is.finite(
      final$loglik
    ) &&
    final_structural$max_abs_gradient <
      1e-4 &&
    final_baseline_score <
      1e-4 &&
    is.finite(
      quadrature_difference
    ) &&
    abs(
      quadrature_difference
    ) <
      0.01

  list(
    fit_ok = fit_ok,
    converged_outer =
      converged,
    converged_final =
      final_converged,
    final_iterations =
      final_iter,
    theta =
      exp(
        q[1]
      ),
    gamma =
      exp(
        q[2]
      ),
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
      final_structural$max_abs_gradient,
    baseline_gradient =
      final_baseline_score,
    selected_method =
      final_structural$selected_method,
    history =
      history,
    weibull_initial_fit =
      weibull_fit
  )
}
