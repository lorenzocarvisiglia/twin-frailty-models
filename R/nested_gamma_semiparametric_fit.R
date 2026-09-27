make_nested_semiparametric_setup <- function(
  sim,
  n_nodes = 64L
) {
  prep <- prepare_nested_fit(
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
    rule = nested_legendre_rule(
      n_nodes
    ),
    n_nodes = n_nodes
  )
}


nested_semiparametric_range_matrix <- function(
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


nested_semiparametric_weights <- function(
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
    nested_semiparametric_range_matrix(
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
    nested_semiparametric_range_matrix(
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
      M =
        M
    )

  list(
    A = A_weight,
    R = A_weight +
      P_weight
  )
}


nested_semiparametric_mz_integral <- function(
  d,
  H,
  theta,
  phi,
  setup,
  reference = FALSE
) {
  d <- as.numeric(d)
  H <- as.numeric(H)

  if (phi <= 1e-6) {
    return(
      nested_log_integral_phi_zero(
        d = d,
        H = H,
        theta = theta
      )
    )
  }

  shape_w <- 1 /
    phi

  constant <-
    d *
      log(phi) +
    lgamma(
      shape_w +
        d
    ) -
    lgamma(
      shape_w
    )

  if (!reference) {
    return(
      nested_log_integral_asymmetric(
        d = d,
        c1 =
          phi *
            H,
        a1 =
          shape_w +
            d,
        c2 =
          rep(
            0,
            length(d)
          ),
        a2 =
          rep(
            0,
            length(d)
          ),
        constant =
          constant,
        theta =
          theta,
        rule =
          setup$rule,
        log_drop =
          20
      )
    )
  }

  vapply(
    seq_along(d),
    function(i) {
      nested_log_integral_reference(
        d = d[i],
        c1 =
          phi *
            H[i],
        a1 =
          shape_w +
            d[i],
        c2 = 0,
        a2 = 0,
        constant =
          constant[i],
        theta =
          theta
      )
    },
    numeric(1)
  )
}


nested_semiparametric_dz_integral <- function(
  d1,
  d2,
  H1,
  H2,
  theta,
  phi,
  setup,
  reference = FALSE
) {
  d1 <- as.numeric(d1)
  d2 <- as.numeric(d2)
  H1 <- as.numeric(H1)
  H2 <- as.numeric(H2)

  d <- d1 +
    d2

  if (phi <= 1e-6) {
    return(
      nested_log_integral_phi_zero(
        d = d,
        H =
          H1 +
          H2,
        theta =
          theta
      )
    )
  }

  shape_w <- 1 /
    phi

  constant <-
    d *
      log(phi) +
    lgamma(
      shape_w +
        d1
    ) +
    lgamma(
      shape_w +
        d2
    ) -
    2 *
      lgamma(
        shape_w
      )

  if (!reference) {
    return(
      nested_log_integral_asymmetric(
        d = d,
        c1 =
          phi *
            H1,
        a1 =
          shape_w +
            d1,
        c2 =
          phi *
            H2,
        a2 =
          shape_w +
            d2,
        constant =
          constant,
        theta =
          theta,
        rule =
          setup$rule,
        log_drop =
          20
      )
    )
  }

  vapply(
    seq_along(d),
    function(i) {
      nested_log_integral_reference(
        d = d[i],
        c1 =
          phi *
            H1[i],
        a1 =
          shape_w +
            d1[i],
        c2 =
          phi *
            H2[i],
        a2 =
          shape_w +
            d2[i],
        constant =
          constant[i],
        theta =
          theta
      )
    },
    numeric(1)
  )
}


nested_semiparametric_components <- function(
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

  theta <- q[1]
  phi <- q[2]
  alpha <- q[3]

  if (
    theta < 0 ||
    phi < 0
  ) {
    return(
      list(
        loglik = -Inf
      )
    )
  }

  weights <-
    nested_semiparametric_weights(
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
      nested_semiparametric_mz_integral(
        d = d,
        H = H_pair,
        theta = theta,
        phi = phi,
        setup = setup,
        reference = reference
      )

    den <-
      nested_semiparametric_mz_integral(
        d =
          rep(
            0,
            length(idx)
          ),
        H = A_pair,
        theta = theta,
        phi = phi,
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
      nested_semiparametric_dz_integral(
        d1 = d1,
        d2 = d2,
        H1 = H[s1],
        H2 = H[s2],
        theta = theta,
        phi = phi,
        setup = setup,
        reference = reference
      )

    den <-
      nested_semiparametric_dz_integral(
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
        theta = theta,
        phi = phi,
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
    H 
    = H,
    weights = weights
  )
}


nested_semiparametric_baseline_denominator <- function(
  setup,
  q,
  delta
) {
  prep <- setup$prep

  theta <- q[1]
  phi <- q[2]

  z <-
    nested_semiparametric_components(
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
      nested_semiparametric_mz_integral(
        d = d,
        H = H_pair,
        theta = theta,
        phi = phi,
        setup = setup
      )

    log_H_next <-
      nested_semiparametric_mz_integral(
        d = d + 1,
        H = H_pair,
        theta = theta,
        phi = phi,
        setup = setup
      )

    log_A <-
      nested_semiparametric_mz_integral(
        d =
          rep(
            0,
            length(idx)
          ),
        H = A_pair,
        theta = theta,
        phi = phi,
        setup = setup
      )

    log_A_next <-
      nested_semiparametric_mz_integral(
        d =
          rep(
            1,
            length(idx)
          ),
        H = A_pair,
        theta = theta,
        phi = phi,
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
      nested_semiparametric_dz_integral(
        d1 = d1,
        d2 = d2,
        H1 = H[s1],
        H2 = H[s2],
        theta = theta,
        phi = phi,
        setup = setup
     )

    log_H1 <-
      nested_semiparametric_dz_integral(
        d1 = d1 + 1,
        d2 = d2,
        H1 = H[s1],
        H2 = H[s2],
        theta = theta,
        phi = phi,
        setup = setup
      )

    log_H2 <-
      nested_semiparametric_dz_integral(
        d1 = d1,
        d2 = d2 + 1,
        H1 = H[s1],
        H2 = H[s2],
        theta = theta,
        phi = phi,
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
      nested_semiparametric_dz_integral(
        d1 = zeros,
        d2 = zeros,
        H1 = A[s1],
        H2 = A[s2],
        theta = theta,
        phi = phi,
        setup = setup
     )

    log_A1 <-
      nested_semiparametric_dz_integral(
        d1 = ones,
        d2 = zeros,
        H1 = A[s1],
        H2 = A[s2],
        theta = theta,
        phi = phi,
        setup = setup
      )

    log_A2 <-
      nested_semiparametric_dz_integral(
        d1 = zeros,
        d2 = ones,
        H1 = A[s1],
        H2 = A[s2],
        theta = theta,
        phi = phi,
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


fit_nested_semiparametric_baseline <- function(
  setup,
  q,
  delta_start,
  max_iter = 200L,
  tolerance = 1e-6
) {
  delta <- as.numeric(
    delta_start
  )

  current <-
    nested_semiparametric_components(
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
      nested_semiparametric_baseline_denominator(
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
        nested_semiparametric_components(
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
      nested_semiparametric_baseline_denominator(
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

fit_nested_gamma_semiparametric_profile <- function(
  sim,
  verbose = FALSE
) {
  setup <-
    make_nested_semiparametric_setup(
      sim
    )

  weibull_fit <-
    fit_nested_gamma(
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
    0.7
  )

  phi0 <- get_est(
    "phi",
    0.7
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
    0,
    -10
  )

  upper <- c(
    10,
    10,
    10
  )

  q0 <- c(
    theta0,
    phi0,
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

    metric <- abs(
      gradient
    )

    lower_hit <-
      q -
        lower <
      1e-6

    upper_hit <-
      upper -
        q <
      1e-6

    if (any(lower_hit)) {
      metric[
        lower_hit
      ] <-
        pmax(
          -gradient[
            lower_hit
          ],
          0
        )
    }

    if (any(upper_hit)) {
      metric[
        upper_hit
      ] <-
        pmax(
          gradient[
            upper_hit
          ],
          0
        )
    }

    max(
      metric
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
        fit_nested_semiparametric_baseline(
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
          fit_nested_semiparametric_baseline(
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
        theta =
          q[1],
        phi =
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
          "| theta",
          round(
            q[1],
            6
          ),
          "| phi",
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
        nested_semiparametric_components(
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

        q_minus[] <-
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
    nested_semiparametric_baseline_denominator(
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
    nested_semiparametric_components(
      setup = setup,
      q = q,
      delta = delta
    )

  reference <-
    nested_semiparametric_components(
      setup = setup,
      q = q,
      delta = delta,
      reference = TRUE
    )

  quadrature_difference <-
    final$loglik -
      reference$loglik

  theta_boundary <-
    q[1] <=
      1e-6

  phi_boundary <-
    q[2] <=
      1e-6

  upper_boundary_hit <-
    any(
      upper -
        q <
        1e-5
    )

  alpha_lower_boundary <-
    q[3] -
      lower[3] <
      1e-5

  boundary_problem <-
    upper_boundary_hit ||
    alpha_lower_boundary

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
    theta =
      q[1],
    phi =
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
    theta_boundary =
      theta_boundary,
    phi_boundary =
      phi_boundary,
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
