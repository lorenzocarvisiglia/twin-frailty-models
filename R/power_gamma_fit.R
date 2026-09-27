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


power_legendre_rule <- function(
    n_nodes
) {

  J <- matrix(
    0,
    nrow = n_nodes,
    ncol = n_nodes
  )

  if (n_nodes > 1L) {

    j <- seq_len(
      n_nodes - 1L
    )

    off <- j /
      sqrt(
        4 * j^2 - 1
      )

    J[cbind(j, j + 1L)] <- off
    J[cbind(j + 1L, j)] <- off
  }

  eig <- eigen(
    J,
    symmetric = TRUE
  )

  ord <- order(
    eig$values
  )

  list(
    nodes = eig$values[
      ord
    ],
    log_weights = log(
      pmax(
        2 *
          eig$vectors[
            1,
            ord
          ]^2,
        .Machine$double.xmin
      )
    )
  )
}


power_mz_mode_logv <- function(
    B,
    d,
    theta,
    gamma
) {

  B <- as.numeric(
    B
  )

  d <- as.numeric(
    d
  )

  k <- 1 /
    theta

  a <- k +
    gamma *
    d

  m <- log(
    theta *
      a
  )

  positive <- B >
    1e-14

  if (any(positive)) {

    root_b <- (
      log(
        a[
          positive
        ]
      ) -
        log(
          B[
            positive
          ] *
            gamma
        )
    ) /
      gamma

    m[
      positive
    ] <- pmin(
      m[
        positive
      ],
      root_b
    )
  }

  for (iter in seq_len(100L)) {

    e1 <- exp(
      pmin(
        m,
        700
      )
    )

    eg <- exp(
      pmin(
        gamma *
          m,
        700
      )
    )

    gradient <- a -
      e1 /
      theta -
      B *
      gamma *
      eg

    curvature <- e1 /
      theta +
      B *
      gamma^2 *
      eg

    curvature <- pmax(
      curvature,
      .Machine$double.xmin
    )

    step <- gradient /
      curvature

    step <- pmax(
      pmin(
        step,
        5
      ),
      -5
    )

    new_m <- m +
      step

    if (
      max(
        abs(
          new_m -
            m
        )
      ) <
        1e-13
    ) {
      m <- new_m
      break
    }

    m <- new_m
  }

  m
}


power_mz_log_integral_reference <- function(
    B,
    d,
    theta,
    gamma
) {

  k <- 1 /
    theta

  if (B <= 1e-14) {
    return(
      gamma *
        d *
        log(theta) +
        lgamma(
          k +
            gamma *
            d
        ) -
        lgamma(k)
    )
  }

  m <- power_mz_mode_logv(
    B = B,
    d = d,
    theta = theta,
    gamma = gamma
  )

  log_kernel <- function(t) {

    (
      k +
        gamma *
        d
    ) *
      t -
      exp(
        pmin(
          t,
          700
        )
      ) /
      theta -
      B *
      exp(
        pmin(
          gamma *
            t,
          700
        )
      ) -
      lgamma(k) -
      k *
      log(theta)
  }

  peak <- log_kernel(
    m
  )

  value <- tryCatch(
    integrate(
      function(t) {

        z <- log_kernel(t) -
          peak

        z[
          !is.finite(z)
        ] <- -Inf

        exp(z)
      },
      lower = m - 50,
      upper = m + 50,
      rel.tol = 1e-11,
      subdivisions = 1000L,
      stop.on.error = TRUE
    )$value,
    error = function(e) {
      NA_real_
    }
  )

  if (
    !is.finite(value) ||
    value <= 0
  ) {
    return(
      NA_real_
    )
  }

  peak +
    log(value)
}


power_mz_log_kernel_vec <- function(
    t,
    B,
    d,
    theta,
    gamma
) {

  k <- 1 /
    theta

  (
    k +
      gamma *
      d
  ) *
    t -
    exp(
      pmin(
        t,
        700
      )
    ) /
    theta -
    B *
    exp(
      pmin(
        gamma *
          t,
        700
      )
    ) -
    lgamma(k) -
    k *
    log(theta)
}


power_mz_find_bounds_vec <- function(
    B,
    d,
    theta,
    gamma,
    log_drop = 20
) {

  m <- power_mz_mode_logv(
    B = B,
    d = d,
    theta = theta,
    gamma = gamma
  )

  peak <- power_mz_log_kernel_vec(
    t = m,
    B = B,
    d = d,
    theta = theta,
    gamma = gamma
  )

  target <- peak -
    log_drop

  left <- m -
    1

  for (iter in seq_len(60L)) {

    g_left <- power_mz_log_kernel_vec(
      t = left,
      B = B,
      d = d,
      theta = theta,
      gamma = gamma
    )

    move <-
      is.finite(g_left) &
      is.finite(target) &
      g_left >
      target

    if (!any(move)) {
      break
    }

    left[
      move
    ] <- m[
      move
    ] -
      2 *
      (
        m[
          move
        ] -
          left[
            move
          ]
      )
  }

  right <- m +
    1

  for (iter in seq_len(60L)) {

    g_right <- power_mz_log_kernel_vec(
      t = right,
      B = B,
      d = d,
      theta = theta,
      gamma = gamma
    )

    move <-
      is.finite(g_right) &
      is.finite(target) &
      g_right >
      target

    if (!any(move)) {
      break
    }

    right[
      move
    ] <- m[
      move
    ] +
      2 *
      (
        right[
          move
        ] -
          m[
            move
          ]
      )
  }

  g_left <- power_mz_log_kernel_vec(
    t = left,
    B = B,
    d = d,
    theta = theta,
    gamma = gamma
  )

  g_right <- power_mz_log_kernel_vec(
    t = right,
    B = B,
    d = d,
    theta = theta,
    gamma = gamma
  )

  valid <-
    is.finite(m) &
    is.finite(target) &
    is.finite(left) &
    is.finite(right) &
    (
      !is.finite(g_left) |
        g_left <=
        target
    ) &
    (
      !is.finite(g_right) |
        g_right <=
        target
    )

  lo <- left
  hi <- m

  for (iter in seq_len(60L)) {

    mid <- (
      lo +
        hi
    ) /
      2

    g_mid <- power_mz_log_kernel_vec(
      t = mid,
      B = B,
      d = d,
      theta = theta,
      gamma = gamma
    )

    below <-
      valid &
      (
        !is.finite(g_mid) |
          g_mid <=
          target
      )

    above <-
      valid &
      !below

    lo[
      below
    ] <- mid[
      below
    ]

    hi[
      above
    ] <- mid[
      above
    ]
  }

  lower <- (
    lo +
      hi
  ) /
    2

  lo <- m
  hi <- right

  for (iter in seq_len(60L)) {

    mid <- (
      lo +
        hi
    ) /
      2

    g_mid <- power_mz_log_kernel_vec(
      t = mid,
      B = B,
      d = d,
      theta = theta,
      gamma = gamma
    )

    above <-
      valid &
      is.finite(g_mid) &
      g_mid >
      target

    below <-
      valid &
      !above

    lo[
      above
    ] <- mid[
      above
    ]

    hi[
      below
    ] <- mid[
      below
    ]
  }

  upper <- (
    lo +
      hi
  ) /
    2

  list(
    lower = lower,
    upper = upper,
    valid = valid
  )
}


power_mz_log_integral_asymmetric <- function(
    B,
    d,
    theta,
    gamma,
    rule,
    log_drop = 20
) {

  B <- as.numeric(
    B
  )

  d <- as.numeric(
    d
  )

  k <- 1 /
    theta

  out <- rep(
    NA_real_,
    length(B)
  )

  zero <- B <=
    1e-14

  if (any(zero)) {

    out[
      zero
    ] <-
      gamma *
      d[
        zero
      ] *
      log(theta) +
      lgamma(
        k +
          gamma *
          d[
            zero
          ]
      ) -
      lgamma(k)
  }

  positive <- !zero

  if (!any(positive)) {
    return(
      out
    )
  }

  Bp <- B[
    positive
  ]

  dp <- d[
    positive
  ]

  bounds <- power_mz_find_bounds_vec(
    B = Bp,
    d = dp,
    theta = theta,
    gamma = gamma,
    log_drop = log_drop
  )

  good <- bounds$valid

  if (!any(good)) {
    return(
      out
    )
  }

  Bg <- Bp[
    good
  ]

  dg <- dp[
    good
  ]

  lower <- bounds$lower[
    good
  ]

  upper <- bounds$upper[
    good
  ]

  midpoint <- (
    lower +
      upper
  ) /
    2

  half_width <- (
    upper -
      lower
  ) /
    2

  tmat <- outer(
    half_width,
    rule$nodes,
    "*"
  )

  tmat <- sweep(
    tmat,
    1,
    midpoint,
    "+"
  )

  a <- k +
    gamma *
    dg

  gmat <- sweep(
    tmat,
    1,
    a,
    "*"
  )

  gmat <-
    gmat -
    exp(
      pmin(
        tmat,
        700
      )
    ) /
    theta

  gmat <-
    gmat -
    sweep(
      exp(
        pmin(
          gamma *
            tmat,
          700
        )
      ),
      1,
      Bg,
      "*"
    )

  gmat <-
    gmat -
    lgamma(k) -
    k *
    log(theta)

  gmat <- sweep(
    gmat,
    2,
    rule$log_weights,
    "+"
  )

  mx <- apply(
    gmat,
    1,
    max
  )

  values <-
    log(
      half_width
    ) +
    mx +
    log(
      rowSums(
        exp(
          sweep(
            gmat,
            1,
            mx,
            "-"
          )
        )
      )
    )

  positive_index <- which(
    positive
  )

  out[
    positive_index[
      good
    ]
  ] <- values

  out
}


make_power_nll <- function(
    prep,
    n_nodes = 64L,
    log_drop = 20,
    reference = FALSE
) {

  rule <- if (!reference) {
    power_legendre_rule(
      n_nodes
    )
  } else {
    NULL
  }

  function(
      q,
      diagnostics = FALSE
  ) {

    q <- as.numeric(
      q
    )

    if (
      length(q) != 5L ||
      any(
        !is.finite(q)
      )
    ) {
      return(
        1e100
      )
    }

    theta <- exp(
      q[1]
    )

    gamma <- exp(
      q[2]
    )

    alpha <- q[3]

    lambda <- exp(
      q[4]
    )

    rho <- exp(
      q[5]
    )

    covariate_factor <- exp(
      alpha *
        prep$y
    )

    interval_exposure <-
      lambda *
      (
        prep$stop^rho -
          prep$start^rho
      ) *
      covariate_factor

    if (
      any(
        !is.finite(
          interval_exposure
        )
      ) ||
      any(
        interval_exposure <
          0
      )
    ) {
      return(
        1e100
      )
    }

    post_pair <- sum_by_group(
      interval_exposure,
      prep$pair_index_long,
      prep$n_pairs
    )

    entry_increment <- pmax(
      prep$entry_interval_end^rho -
        prep$entry_interval_start^rho,
      0
    )

    entry_subject <-
      lambda *
      rowSums(
        entry_increment *
          exp(
            alpha *
              prep$y_history
          )
      )

    if (
      any(
        !is.finite(
          entry_subject
        )
      ) ||
      any(
        entry_subject <
          0
      )
    ) {
      return(
        1e100
      )
    }

    entry_pair <- sum_by_group(
      entry_subject,
      prep$pair_index_subject,
      prep$n_pairs
    )

    total_pair <-
      entry_pair +
      post_pair

    event_index <-
      prep$status ==
      1L

    log_hazard <- 0

    if (any(event_index)) {

      log_hazard <- sum(
        log(lambda) +
          log(rho) +
          (
            rho -
              1
          ) *
          prep$log_stop[
            event_index
          ] +
          alpha *
          prep$y[
            event_index
          ]
      )
    }

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

      A <- entry_pair[
        idx
      ]

      H <- total_pair[
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
              A
          ) -
          (
            shape +
              d
          ) *
          log1p(
            theta *
              H
          )
      )
    }

    loglik_mz <- 0
    mz_integral_count <- 0L

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

      A <- entry_pair[
        idx
      ]

      H <- total_pair[
        idx
      ]

      if (reference) {

        num <- vapply(
          seq_along(H),
          function(i) {
            power_mz_log_integral_reference(
              B = H[i],
              d = d[i],
              theta = theta,
              gamma = gamma
            )
          },
          numeric(1)
        )

        den <- vapply(
          seq_along(A),
          function(i) {
            power_mz_log_integral_reference(
              B = A[i],
              d = 0,
              theta = theta,
              gamma = gamma
            )
          },
          numeric(1)
        )

      } else {

        num <- power_mz_log_integral_asymmetric(
          B = H,
          d = d,
          theta = theta,
          gamma = gamma,
          rule = rule,
          log_drop = log_drop
        )

        den <- power_mz_log_integral_asymmetric(
          B = A,
          d = rep(
            0,
            length(A)
          ),
          theta = theta,
          gamma = gamma,
          rule = rule,
          log_drop = log_drop
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
          1e100
        )
      }

      loglik_mz <- sum(
        num -
          den
      )

      mz_integral_count <-
        length(num) +
        length(den)
    }

    loglik <-
      log_hazard +
      loglik_dz +
      loglik_mz

    if (
      !is.finite(
        loglik
      )
    ) {
      return(
        1e100
      )
    }

    nll_value <- -loglik

    if (!diagnostics) {
      return(
        nll_value
      )
    }

    list(
      nll = nll_value,
      loglik = loglik,
      loglik_mz = loglik_mz,
      mz_integral_count =
        mz_integral_count
    )
  }
}


fit_power_gamma <- function(
    sim
) {

  prep <- prepare_power_fit(
    sim
  )

  nll <- make_power_nll(
    prep,
    n_nodes = 64L,
    log_drop = 20,
    reference = FALSE
  )

  lower <- c(
    log(0.02),
    log(0.05),
    -10,
    log(1e-5),
    log(0.3)
  )

  upper <- c(
    log(10),
    log(8),
    10,
    log(2),
    log(8)
  )

  start1 <- c(
    log(0.5),
    log(1),
    0,
    log(0.1),
    log(2)
  )

  fit1 <- tryCatch(
    nlminb(
      start = start1,
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
    !is.null(fit1) &&
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

    fit2 <- tryCatch(
      optim(
        par = fit1$par,
        fn = nll,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(
          maxit = 1000,
          factr = 1e7,
          pgtol = 1e-8
        )
      ),
      error = function(e) {
        NULL
      }
    )

    if (
      !is.null(fit2) &&
      is.finite(
        fit2$value
      )
    ) {

      candidates$L_BFGS_B <- list(
        par = fit2$par,
        objective = fit2$value,
        convergence = fit2$convergence,
        message = fit2$message
      )
    }
  }

  has_good <- if (
    length(candidates) >
      0L
  ) {

    any(
      vapply(
        candidates,
        function(x) {
          x$convergence ==
            0L
        },
        logical(1)
      )
    )

  } else {

    FALSE
  }

  if (!has_good) {

    start2 <- c(
      log(1.2),
      log(2),
      -1,
      log(0.05),
      log(2.5)
    )

    fit3 <- tryCatch(
      nlminb(
        start = start2,
        objective = nll,
        lower = lower,
        upper = upper,
        control = list(
          eval.max = 1800,
          iter.max = 1200,
          rel.tol = 1e-9,
          x.tol = 1e-9
        )
      ),
      error = function(e) {
        NULL
      }
    )

    if (
      !is.null(fit3) &&
      is.finite(
        fit3$objective
      )
    ) {

      candidates$fallback <- list(
        par = fit3$par,
        objective = fit3$objective,
        convergence = fit3$convergence,
        message = fit3$message
      )
    }
  }

  if (length(candidates) == 0L) {
    stop(
      "all power gamma optimizers failed"
    )
  }

  convergence_zero <- vapply(
    candidates,
    function(x) {
      x$convergence ==
        0L
    },
    logical(1)
  )

  candidate_objectives <- vapply(
    candidates,
    function(x) {
      x$objective
    },
    numeric(1)
  )

  if (any(convergence_zero)) {

    good_idx <- which(
      convergence_zero
    )

    selected_idx <- good_idx[
      which.min(
        candidate_objectives[
          good_idx
        ]
      )
    ]

  } else {

    selected_idx <- which.min(
      candidate_objectives
    )
  }

  fit <- candidates[[
    selected_idx
  ]]

  selected_method <- names(
    candidates
  )[
    selected_idx
  ]

  q_hat <- fit$par

  pre_refinement_gradient <- tryCatch(
    bounded_gradient(
      fn = nll,
      x = q_hat,
      lower = lower,
      upper = upper,
      rel_step = 1e-4
    ),
    error = function(e) {
      rep(
        NA_real_,
        5
      )
    }
  )

  pre_refinement_max_abs_gradient <- if (
    all(
      is.finite(
        pre_refinement_gradient
      )
    )
  ) {

    max(
      abs(
        pre_refinement_gradient
      )
    )

  } else {

    NA_real_
  }

  strict_refinement_attempted <-
    !is.finite(
      pre_refinement_max_abs_gradient
    ) ||
    pre_refinement_max_abs_gradient >
    0.005

  if (strict_refinement_attempted) {

    refinement_candidates <- list(
      current = fit
    )

    fit_strict_nl <- tryCatch(
      nlminb(
        start = q_hat,
        objective = nll,
        lower = lower,
        upper = upper,
        control = list(
          eval.max = 4000,
          iter.max = 2500,
          rel.tol = 1e-12,
          x.tol = 1e-12
        )
      ),
      error = function(e) {
        NULL
      }
    )

    if (
      !is.null(
        fit_strict_nl
      ) &&
      is.finite(
        fit_strict_nl$objective
      )
    ) {

      refinement_candidates$strict_nlminb <- list(
        par = fit_strict_nl$par,
        objective = fit_strict_nl$objective,
        convergence = fit_strict_nl$convergence,
        message = fit_strict_nl$message
      )
    }

    strict_start <- if (
      !is.null(
        fit_strict_nl
      ) &&
      is.finite(
        fit_strict_nl$objective
      )
    ) {

      fit_strict_nl$par

    } else {

      q_hat
    }

    fit_strict_lb <- tryCatch(
      optim(
        par = strict_start,
        fn = nll,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(
          maxit = 3000,
          factr = 1e4,
          pgtol = 1e-10
        )
      ),
      error = function(e) {
        NULL
      }
    )

    if (
      !is.null(
        fit_strict_lb
      ) &&
      is.finite(
        fit_strict_lb$value
      )
    ) {

      refinement_candidates$strict_L_BFGS_B <- list(
        par = fit_strict_lb$par,
        objective = fit_strict_lb$value,
        convergence = fit_strict_lb$convergence,
        message = fit_strict_lb$message
      )
    }

    refinement_convergence_zero <- vapply(
      refinement_candidates,
      function(x) {
        x$convergence ==
          0L
      },
      logical(1)
    )

    refinement_objectives <- vapply(
      refinement_candidates,
      function(x) {
        x$objective
      },
      numeric(1)
    )

    if (any(refinement_convergence_zero)) {

      good_idx <- which(
        refinement_convergence_zero
      )

      refinement_selected_idx <- good_idx[
        which.min(
          refinement_objectives[
            good_idx
          ]
        )
      ]

    } else {

      refinement_selected_idx <- which.min(
        refinement_objectives
      )
    }

    fit <- refinement_candidates[[
      refinement_selected_idx
    ]]

    selected_method <- names(
      refinement_candidates
    )[
      refinement_selected_idx
    ]

    q_hat <- fit$par
  }

  estimate <- c(
    theta = exp(
      q_hat[1]
    ),
    gamma = exp(
      q_hat[2]
    ),
    alpha = q_hat[3],
    lambda = exp(
      q_hat[4]
    ),
    rho = exp(
      q_hat[5]
    )
  )

  distance_lower <- q_hat -
    lower

  distance_upper <- upper -
    q_hat

  boundary_hit <- any(
    distance_lower <
      1e-5 |
      distance_upper <
      1e-5
  )

  gradient <- tryCatch(
    bounded_gradient(
      fn = nll,
      x = q_hat,
      lower = lower,
      upper = upper,
      rel_step = 1e-4
    ),
    error = function(e) {
      rep(
        NA_real_,
        5
      )
    }
  )

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

    NA_real_
  }

  asymmetric_details <- nll(
    q_hat,
    diagnostics = TRUE
  )

  nll_reference <- make_power_nll(
    prep,
    reference = TRUE
  )

  reference_details <- nll_reference(
    q_hat,
    diagnostics = TRUE
  )

  loglik_asymmetric <-
    asymmetric_details$loglik

  loglik_reference <-
    reference_details$loglik

  quadrature_difference <-
    loglik_asymmetric -
    loglik_reference

  fit_ok <-
    fit$convergence ==
    0L &&
    !boundary_hit &&
    is.finite(
      fit$objective
    ) &&
    is.finite(
      max_abs_gradient
    ) &&
    max_abs_gradient <
    0.01 &&
    is.finite(
      quadrature_difference
    ) &&
    abs(
      quadrature_difference
    ) <
    0.01

  list(
    estimates = estimate,
    convergence = fit$convergence,
    fit_ok = fit_ok,
    message = paste(
      selected_method,
      ifelse(
        is.null(
          fit$message
        ),
        "",
        as.character(
          fit$message
        )
      )
    ),
    logLik = loglik_asymmetric,
    n_events = sim$n_events,
    event_proportion =
      sim$event_proportion,
    diagnostics = list(
      selected_method =
        selected_method,
      objective =
        fit$objective,
      boundary_hit =
        boundary_hit,
      pre_refinement_max_abs_gradient =
        pre_refinement_max_abs_gradient,
      strict_refinement_attempted =
        strict_refinement_attempted,
      max_abs_gradient =
        max_abs_gradient,
      asymmetric_logLik =
        loglik_asymmetric,
      reference_logLik =
        loglik_reference,
      quadrature_difference =
        quadrature_difference,
      quadrature_nodes = 64L,
      quadrature_log_drop = 20,
      mz_integral_count =
        asymmetric_details$mz_integral_count,
      n_pairs_generated =
        sim$n_pairs_generated,
      n_pairs_observed =
        sim$n_pairs_observed,
      n_dz_observed =
        sim$n_dz_observed,
      n_mz_observed =
        sim$n_mz_observed,
      n_subjects_observed =
        sim$n_subjects_observed,
      generated_frailty_mean =
        sim$generated_frailty_mean,
      generated_frailty_variance =
        sim$generated_frailty_variance
    )
  )
}
