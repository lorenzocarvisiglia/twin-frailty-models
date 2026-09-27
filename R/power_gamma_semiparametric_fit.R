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

fit_power_gamma_semiparametric_profile <- function(
  sim,
  verbose = FALSE
) {
  setup <- make_power_semiparametric_setup(sim)
  weibull_fit <- fit_power_gamma(sim = sim)

  get_est <- function(name, fallback) {
    z <- weibull_fit$estimates[[name]]
    if (length(z) == 1L && is.finite(z)) {
      as.numeric(z)
    } else {
      fallback
    }
  }

  theta0 <- get_est("theta", 0.5)
  gamma0 <- get_est("gamma", 1)
  alpha0 <- get_est("alpha", 0)
  lambda0 <- get_est("lambda", 0.1)
  rho0 <- get_est("rho", 2)

  q0 <- c(log(theta0), log(gamma0), alpha0)
  previous_times <- c(0, setup$event_times[-setup$M])
  delta0 <- lambda0 * (
    setup$event_times^rho0 -
      previous_times^rho0
  )
  delta0 <- pmax(delta0, 1e-10)

  lower <- c(log(0.02), log(0.05), -10)
  upper <- c(log(10), log(8), 10)

  make_solver <- function(
    delta_initial,
    baseline_tolerance,
    baseline_max_iter,
    phase
  ) {
    state <- new.env(parent = emptyenv())
    state$q <- NULL
    state$delta <- delta_initial
    state$profile <- NULL
    state$n_profile <- 0L
    state$history <- list()

    same_q <- function(x, y) {
      !is.null(y) &&
        length(x) == length(y) &&
        max(abs(x - y)) < 1e-13
    }

    profile_at <- function(q) {
      q <- as.numeric(q)

      if (same_q(q, state$q)) {
        return(state$profile)
      }

      baseline_fit <- tryCatch(
        fit_power_semiparametric_baseline(
          setup = setup,
          q = q,
          delta_start = state$delta,
          max_iter = baseline_max_iter,
          tolerance = baseline_tolerance
        ),
        error = function(e) NULL
      )

      if (
        is.null(baseline_fit) ||
        !is.finite(baseline_fit$loglik)
      ) {
        baseline_fit <- tryCatch(
          fit_power_semiparametric_baseline(
            setup = setup,
            q = q,
            delta_start = delta_initial,
            max_iter = baseline_max_iter,
            tolerance = baseline_tolerance
          ),
          error = function(e) NULL
        )
      }

      if (
        is.null(baseline_fit) ||
        !is.finite(baseline_fit$loglik)
      ) {
        return(NULL)
      }

      state$n_profile <- state$n_profile + 1L

      out <- list(
        q = q,
        delta = baseline_fit$delta,
        loglik = baseline_fit$loglik,
        nll = -baseline_fit$loglik,
        baseline_iterations = baseline_fit$iterations,
        baseline_gradient = baseline_fit$max_eta_score,
        baseline_convergence = baseline_fit$convergence
      )

      state$q <- q
      state$delta <- baseline_fit$delta
      state$profile <- out
      state$history[[length(state$history) + 1L]] <- list(
        phase = phase,
        evaluation = state$n_profile,
        theta = exp(q[1]),
        gamma = exp(q[2]),
        alpha = q[3],
        logLik = baseline_fit$loglik,
        baseline_iterations = baseline_fit$iterations,
        baseline_gradient = baseline_fit$max_eta_score
      )

      if (verbose) {
        cat(
          phase,
          state$n_profile,
          "| theta",
          round(exp(q[1]), 6),
          "| gamma",
          round(exp(q[2]), 6),
          "| alpha",
          round(q[3], 6),
          "| logLik",
          round(baseline_fit$loglik, 6),
          "| baseline iter",
          baseline_fit$iterations,
          "| baseline score",
          signif(baseline_fit$max_eta_score, 4),
          "\n"
        )
      }

      out
    }

    fixed_delta_nll <- function(q, delta) {
      z <- power_semiparametric_components(
        setup = setup,
        q = q,
        delta = delta
      )

      if (!is.finite(z$loglik)) {
        return(1e100)
      }

      -z$loglik
    }

    envelope_gradient <- function(
      q,
      delta,
      rel_step = 1e-4
    ) {
      g <- numeric(length(q))

      for (k in seq_along(q)) {
        h <- rel_step * max(1, abs(q[k]))
        q_minus <- q
        q_plus <- q

        q_minus[k] <- max(
          lower[k],
          q[k] - h
        )

        q_plus[k] <- min(
          upper[k],
          q[k] + h
        )

        if (q_plus[k] > q_minus[k]) {
          f_minus <- fixed_delta_nll(
            q_minus,
            delta
          )
          f_plus <- fixed_delta_nll(
            q_plus,
            delta
          )

          g[k] <- (
            f_plus -
              f_minus
          ) / (
            q_plus[k] -
              q_minus[k]
          )
        } else {
          g[k] <- 0
        }
      }

      g
    }

    list(
      profile_at = profile_at,
      envelope_gradient = envelope_gradient,
      state = state
    )
  }

  run_optimizer <- function(
    solver,
    q_start,
    factr,
    pgtol
  ) {
    initial <- solver$profile_at(q_start)

    if (is.null(initial)) {
      stop("initial profile failed")
    }

    reference_nll <- initial$nll

    objective <- function(q) {
      z <- solver$profile_at(q)

      if (is.null(z)) {
        anchor <- if (is.null(solver$state$q)) {
          q_start
        } else {
          solver$state$q
        }

        dq <- q - anchor

        return(
          1e6 +
            1e4 *
              sum(dq^2)
        )
      }

      z$nll -
        reference_nll
    }

    gradient <- function(q) {
      z <- solver$profile_at(q)

      if (is.null(z)) {
        anchor <- if (is.null(solver$state$q)) {
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

  main_solver <- make_solver(
    delta_initial = delta0,
    baseline_tolerance = 1e-7,
    baseline_max_iter = 1000L,
    phase = "profile"
  )

  main_fit <- run_optimizer(
    solver = main_solver,
    q_start = q0,
    factr = 100,
    pgtol = 1e-7
  )

  main_profile <- main_solver$profile_at(
    main_fit$par
  )

  if (is.null(main_profile)) {
    stop("final main profile failed")
  }

  main_gradient <- main_solver$envelope_gradient(
    q = main_fit$par,
    delta = main_profile$delta,
    rel_step = 1e-5
  )

  main_max_gradient <- max(
    abs(main_gradient)
  )

  polished <- main_max_gradient > 1e-4

  if (polished) {
    polish_solver <- make_solver(
      delta_initial = main_profile$delta,
      baseline_tolerance = 1e-9,
      baseline_max_iter = 2000L,
      phase = "polish"
    )

    polish_fit <- run_optimizer(
      solver = polish_solver,
      q_start = main_fit$par,
      factr = 1,
      pgtol = 1e-9
    )

    final_fit <- polish_fit
    final_solver <- polish_solver
    final_profile <- polish_solver$profile_at(
      polish_fit$par
    )

    history <- c(
      main_solver$state$history,
      polish_solver$state$history
    )

    profile_evaluations <-
      main_solver$state$n_profile +
      polish_solver$state$n_profile

    polishing_evaluations <-
      polish_solver$state$n_profile
  } else {
    final_fit <- main_fit
    final_solver <- main_solver
    final_profile <- main_profile
    history <- main_solver$state$history
    profile_evaluations <- main_solver$state$n_profile
    polishing_evaluations <- 0L
  }

  if (is.null(final_profile)) {
    stop("final profile failed")
  }

  q <- final_fit$par
  delta <- final_profile$delta

  final_gradient <- final_solver$envelope_gradient(
    q = q,
    delta = delta,
    rel_step = 1e-5
  )

  structural_gradient <- max(
    abs(final_gradient)
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

  final <- power_semiparametric_components(
    setup = setup,
    q = q,
    delta = delta
  )

  reference <- power_semiparametric_components(
    setup = setup,
    q = q,
    delta = delta,
    reference = TRUE
  )

  quadrature_difference <-
    final$loglik -
      reference$loglik

  boundary_hit <- any(
    q -
      lower <
      1e-5 |
      upper -
        q <
        1e-5
  )

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
    !boundary_hit &&
    is.finite(
      quadrature_difference
    ) &&
    abs(
      quadrature_difference
    ) <
      0.01

  list(
    fit_ok = fit_ok,
    converged_outer = converged_outer,
    converged_final = converged_final,
    final_iterations = profile_evaluations,
    theta = exp(q[1]),
    gamma = exp(q[2]),
    alpha = q[3],
    delta = delta,
    event_times = setup$event_times,
    logLik = final$loglik,
    reference_logLik = reference$loglik,
    quadrature_difference = quadrature_difference,
    structural_gradient = structural_gradient,
    baseline_gradient = final_baseline_score,
    selected_method = if (polished) {
      "profile_L_BFGS_B_polished"
    } else {
      "profile_L_BFGS_B"
    },
    history = history,
    weibull_initial_fit = weibull_fit,
    profile_evaluations = profile_evaluations,
    polishing_evaluations = polishing_evaluations,
    polished = polished,
    main_structural_gradient = main_max_gradient,
    optimizer_convergence = final_fit$convergence,
    optimizer_message = final_fit$message
  )
}

fit_power_gamma_semiparametric <- function(
  sim,
  verbose = FALSE
) {
  fit_power_gamma_semiparametric_profile(
    sim = sim,
    verbose = verbose
  )
}
