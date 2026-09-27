prepare_nested_fit <- function(
    sim
) {

  prep <- prepare_frailty_fit(
    sim
  )

  subject1 <- seq(
    1L,
    prep$n_subjects,
    by = 2L
  )

  subject2 <- subject1 +
    1L

  if (
    !all(
      prep$subjects$pair_index[
        subject1
      ] ==
        prep$subjects$pair_index[
          subject2
        ]
    )
  ) {
    stop(
      "subjects are not correctly paired"
    )
  }

  prep$subject1 <- subject1
  prep$subject2 <- subject2

  prep$subject_status <-
    prep$subjects$status

  prep
}


nested_legendre_rule <- function(
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
        4 * j^2 -
          1
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


nested_log_kernel_vec <- function(
    t,
    d,
    c1,
    a1,
    c2,
    a2,
    constant,
    theta
) {

  k <- 1 /
    theta

  ev <- exp(
    pmin(
      t,
      700
    )
  )

  (
    k +
      d
  ) *
    t -
    ev /
    theta -
    a1 *
    log1p(
      c1 *
        ev
    ) -
    a2 *
    log1p(
      c2 *
        ev
    ) -
    lgamma(k) -
    k *
    log(theta) +
    constant
}


nested_mode_vec <- function(
    d,
    c1,
    a1,
    c2,
    a2,
    theta
) {

  k <- 1 /
    theta

  m <- log(
    theta *
      (
        k +
          d
      )
  )

  for (iter in seq_len(100L)) {

    ev <- exp(
      pmin(
        m,
        700
      )
    )

    x1 <- c1 *
      ev

    x2 <- c2 *
      ev

    p1 <- x1 /
      (
        1 +
          x1
      )

    p2 <- x2 /
      (
        1 +
          x2
      )

    p1[
      is.infinite(x1)
    ] <- 1

    p2[
      is.infinite(x2)
    ] <- 1

    gradient <-
      k +
      d -
      ev /
      theta -
      a1 *
      p1 -
      a2 *
      p2

    curvature <-
      ev /
      theta +
      a1 *
      p1 *
      (
        1 -
          p1
      ) +
      a2 *
      p2 *
      (
        1 -
          p2
      )

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


nested_find_bounds_vec <- function(
    d,
    c1,
    a1,
    c2,
    a2,
    constant,
    theta,
    log_drop = 20
) {

  m <- nested_mode_vec(
    d = d,
    c1 = c1,
    a1 = a1,
    c2 = c2,
    a2 = a2,
    theta = theta
  )

  peak <- nested_log_kernel_vec(
    t = m,
    d = d,
    c1 = c1,
    a1 = a1,
    c2 = c2,
    a2 = a2,
    constant = constant,
    theta = theta
  )

  target <- peak -
    log_drop

  left <- m -
    1

  for (iter in seq_len(60L)) {

    g_left <- nested_log_kernel_vec(
      t = left,
      d = d,
      c1 = c1,
      a1 = a1,
      c2 = c2,
      a2 = a2,
      constant = constant,
      theta = theta
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

    g_right <- nested_log_kernel_vec(
      t = right,
      d = d,
      c1 = c1,
      a1 = a1,
      c2 = c2,
      a2 = a2,
      constant = constant,
      theta = theta
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

  g_left <- nested_log_kernel_vec(
    t = left,
    d = d,
    c1 = c1,
    a1 = a1,
    c2 = c2,
    a2 = a2,
    constant = constant,
    theta = theta
  )

  g_right <- nested_log_kernel_vec(
    t = right,
    d = d,
    c1 = c1,
    a1 = a1,
    c2 = c2,
    a2 = a2,
    constant = constant,
    theta = theta
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

    g_mid <- nested_log_kernel_vec(
      t = mid,
      d = d,
      c1 = c1,
      a1 = a1,
      c2 = c2,
      a2 = a2,
      constant = constant,
      theta = theta
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

    g_mid <- nested_log_kernel_vec(
      t = mid,
      d = d,
      c1 = c1,
      a1 = a1,
      c2 = c2,
      a2 = a2,
      constant = constant,
      theta = theta
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


nested_log_integral_asymmetric <- function(
    d,
    c1,
    a1,
    c2,
    a2,
    constant,
    theta,
    rule,
    log_drop = 20
) {

  d <- as.numeric(d)
  c1 <- as.numeric(c1)
  a1 <- as.numeric(a1)
  c2 <- as.numeric(c2)
  a2 <- as.numeric(a2)
  constant <- as.numeric(constant)

  if (theta <= 1e-6) {
    return(
      constant -
        a1 *
        log1p(c1) -
        a2 *
        log1p(c2)
    )
  }

  out <- rep(
    NA_real_,
    length(d)
  )

  bounds <- nested_find_bounds_vec(
    d = d,
    c1 = c1,
    a1 = a1,
    c2 = c2,
    a2 = a2,
    constant = constant,
    theta = theta,
    log_drop = log_drop
  )

  good <- bounds$valid

  if (!any(good)) {
    return(
      out
    )
  }

  dg <- d[
    good
  ]

  c1g <- c1[
    good
  ]

  a1g <- a1[
    good
  ]

  c2g <- c2[
    good
  ]

  a2g <- a2[
    good
  ]

  constantg <- constant[
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

  ev <- exp(
    pmin(
      tmat,
      700
    )
  )

  k <- 1 /
    theta

  gmat <- sweep(
    tmat,
    1,
    k +
      dg,
    "*"
  )

  gmat <- gmat -
    ev /
    theta

  term1 <- log1p(
    sweep(
      ev,
      1,
      c1g,
      "*"
    )
  )

  term2 <- log1p(
    sweep(
      ev,
      1,
      c2g,
      "*"
    )
  )

  term1 <- sweep(
    term1,
    1,
    a1g,
    "*"
  )

  term2 <- sweep(
    term2,
    1,
    a2g,
    "*"
  )

  gmat <-
    gmat -
    term1 -
    term2 -
    lgamma(k) -
    k *
    log(theta)

  gmat <- sweep(
    gmat,
    1,
    constantg,
    "+"
  )

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

  out[
    good
  ] <- values

  out
}


nested_log_integral_reference <- function(
    d,
    c1,
    a1,
    c2,
    a2,
    constant,
    theta
) {

  if (theta <= 1e-6) {
    return(
      constant -
        a1 *
        log1p(c1) -
        a2 *
        log1p(c2)
    )
  }

  m <- nested_mode_vec(
    d = d,
    c1 = c1,
    a1 = a1,
    c2 = c2,
    a2 = a2,
    theta = theta
  )

  peak <- nested_log_kernel_vec(
    t = m,
    d = d,
    c1 = c1,
    a1 = a1,
    c2 = c2,
    a2 = a2,
    constant = constant,
    theta = theta
  )

  value <- tryCatch(
    integrate(
      function(t) {

        z <- nested_log_kernel_vec(
          t = t,
          d = d,
          c1 = c1,
          a1 = a1,
          c2 = c2,
          a2 = a2,
          constant = constant,
          theta = theta
        ) -
          peak

        z[
          !is.finite(z)
        ] <- -Inf

        exp(z)
      },
      lower = -Inf,
      upper = Inf,
      rel.tol = 1e-10,
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


nested_log_integral_phi_zero <- function(
    d,
    H,
    theta
) {

  d <- as.numeric(
    d
  )

  H <- as.numeric(
    H
  )

  if (theta <= 1e-6) {
    return(
      -H
    )
  }

  k <- 1 /
    theta

  d *
    log(theta) +
    lgamma(
      k +
        d
    ) -
    lgamma(k) -
    (
      k +
        d
    ) *
    log1p(
      theta *
        H
    )
}


make_nested_nll <- function(
    prep,
    n_nodes = 64L,
    log_drop = 20,
    reference = FALSE
) {

  rule <- if (!reference) {
    nested_legendre_rule(
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

    theta <- q[1]
    phi <- q[2]
    alpha <- q[3]

    lambda <- exp(
      q[4]
    )

    rho <- exp(
      q[5]
    )

    if (
      theta < 0 ||
      phi < 0 ||
      !is.finite(phi) ||
      !is.finite(lambda) ||
      !is.finite(rho)
    ) {
      return(
        1e100
      )
    }

    interval_exposure <-
      lambda *
      (
        prep$stop^rho -
          prep$start^rho
      ) *
      exp(
        alpha *
          prep$y
      )

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

    post_subject <- sum_by_group(
      interval_exposure,
      prep$long_subject_index,
      prep$n_subjects
    )

    entry_subject <-
      lambda *
      prep$entry^rho *
      exp(
        alpha *
          prep$y_entry
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

    total_subject <-
      entry_subject +
      post_subject

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

    phi_boundary_eval <-
      phi <=
      1e-6

    shape_w <- if (
      phi_boundary_eval
    ) {
      NA_real_
    } else {
      1 /
        phi
    }

    loglik_mz <- 0
    loglik_dz <- 0
    v_integral_count <- 0L

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

      A <-
        entry_subject[
          s1
        ] +
        entry_subject[
          s2
        ]

      H <-
        total_subject[
          s1
        ] +
        total_subject[
          s2
        ]

      if (phi_boundary_eval) {

        num <- nested_log_integral_phi_zero(
          d = d,
          H = H,
          theta = theta
        )

        den <- nested_log_integral_phi_zero(
          d = rep(
            0,
            length(d)
          ),
          H = A,
          theta = theta
        )

      } else {

        constant_num <-
          d *
          log(phi) +
          lgamma(
            shape_w +
              d
          ) -
          lgamma(
            shape_w
          )

        if (reference) {

          num <- vapply(
            seq_along(d),
            function(i) {
              nested_log_integral_reference(
                d = d[i],
                c1 = phi *
                  H[i],
                a1 = shape_w +
                  d[i],
                c2 = 0,
                a2 = 0,
                constant =
                  constant_num[i],
                theta = theta
              )
            },
            numeric(1)
          )

          den <- vapply(
            seq_along(d),
            function(i) {
              nested_log_integral_reference(
                d = 0,
                c1 = phi *
                  A[i],
                a1 = shape_w,
                c2 = 0,
                a2 = 0,
                constant = 0,
                theta = theta
              )
            },
            numeric(1)
          )

        } else {

          num <- nested_log_integral_asymmetric(
            d = d,
            c1 = phi *
              H,
            a1 = shape_w +
              d,
            c2 = rep(
              0,
              length(d)
            ),
            a2 = rep(
              0,
              length(d)
            ),
            constant =
              constant_num,
            theta = theta,
            rule = rule,
            log_drop = log_drop
          )

          den <- nested_log_integral_asymmetric(
            d = rep(
              0,
              length(d)
            ),
            c1 = phi *
              A,
            a1 = rep(
              shape_w,
              length(d)
            ),
            c2 = rep(
              0,
              length(d)
            ),
            a2 = rep(
              0,
              length(d)
            ),
            constant = rep(
              0,
              length(d)
            ),
            theta = theta,
            rule = rule,
            log_drop = log_drop
          )
        }
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

      v_integral_count <-
        v_integral_count +
        length(num) +
        length(den)
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

      d <- d1 +
        d2

      A1 <- entry_subject[
        s1
      ]

      A2 <- entry_subject[
        s2
      ]

      H1 <- total_subject[
        s1
      ]

      H2 <- total_subject[
        s2
      ]

      if (phi_boundary_eval) {

        num <- nested_log_integral_phi_zero(
          d = d,
          H =
            H1 +
            H2,
          theta = theta
        )

        den <- nested_log_integral_phi_zero(
          d = rep(
            0,
            length(d)
          ),
          H =
            A1 +
            A2,
          theta = theta
        )

      } else {

        constant_num <-
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

        if (reference) {

          num <- vapply(
            seq_along(d),
            function(i) {
              nested_log_integral_reference(
                d = d[i],
                c1 = phi *
                  H1[i],
                a1 = shape_w +
                  d1[i],
                c2 = phi *
                  H2[i],
                a2 = shape_w +
                  d2[i],
                constant =
                  constant_num[i],
                theta = theta
              )
            },
            numeric(1)
          )

          den <- vapply(
            seq_along(d),
            function(i) {
              nested_log_integral_reference(
                d = 0,
                c1 = phi *
                  A1[i],
                a1 = shape_w,
                c2 = phi *
                  A2[i],
                a2 = shape_w,
                constant = 0,
                theta = theta
              )
            },
            numeric(1)
          )

        } else {

          num <- nested_log_integral_asymmetric(
            d = d,
            c1 = phi *
              H1,
            a1 = shape_w +
              d1,
            c2 = phi *
              H2,
            a2 = shape_w +
              d2,
            constant =
              constant_num,
            theta = theta,
            rule = rule,
            log_drop = log_drop
          )

          den <- nested_log_integral_asymmetric(
            d = rep(
              0,
              length(d)
            ),
            c1 = phi *
              A1,
            a1 = rep(
              shape_w,
              length(d)
            ),
            c2 = phi *
              A2,
            a2 = rep(
              shape_w,
              length(d)
            ),
            constant = rep(
              0,
              length(d)
            ),
            theta = theta,
            rule = rule,
            log_drop = log_drop
          )
        }
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

      loglik_dz <- sum(
        num -
          den
      )

      v_integral_count <-
        v_integral_count +
        length(num) +
        length(den)
    }

    loglik <-
      log_hazard +
      loglik_mz +
      loglik_dz

    if (
      !is.finite(
        loglik
      )
    ) {
      return(
        1e100
      )
    }

    nll_value <-
      -loglik

    if (!diagnostics) {
      return(
        nll_value
      )
    }

    list(
      nll = nll_value,
      loglik = loglik,
      loglik_mz = loglik_mz,
      loglik_dz = loglik_dz,
      v_integral_count =
        v_integral_count
    )
  }
}


fit_nested_gamma <- function(
    sim
) {

  prep <- prepare_nested_fit(
    sim
  )

  nll <- make_nested_nll(
    prep,
    n_nodes = 64L,
    log_drop = 20,
    reference = FALSE
  )

  lower <- c(
    0,
    0,
    -10,
    log(1e-5),
    log(0.3)
  )

  upper <- c(
    10,
    10,
    10,
    log(2),
    log(8)
  )

  start1 <- c(
    0.7,
    0.7,
    0,
    log(0.05),
    log(2)
  )

  fit1 <- tryCatch(
    nlminb(
      start = start1,
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
          maxit = 1200,
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
      1.2,
      0.5,
      -1,
      log(0.02),
      log(2.5)
    )

    fit3 <- tryCatch(
      nlminb(
        start = start2,
        objective = nll,
        lower = lower,
        upper = upper,
        control = list(
          eval.max = 2200,
          iter.max = 1500,
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
      "all nested gamma optimizers failed"
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

  pre_theta_boundary <-
    q_hat[1] <=
    1e-6

  pre_phi_boundary <-
    q_hat[2] <=
    1e-6

  pre_refinement_gradient_metric <- if (
    all(
      is.finite(
        pre_refinement_gradient
      )
    )
  ) {

    metric_components <- abs(
      pre_refinement_gradient
    )

    if (pre_theta_boundary) {
      metric_components[1] <- max(
        -pre_refinement_gradient[1],
        0
      )
    }

    if (pre_phi_boundary) {
      metric_components[2] <- max(
        -pre_refinement_gradient[2],
        0
      )
    }

    max(
      metric_components
    )

  } else {

    NA_real_
  }

  strict_refinement_attempted <-
    !is.finite(
      pre_refinement_gradient_metric
    ) ||
    pre_refinement_gradient_metric >
    0.005

  if (strict_refinement_attempted) {

    strict_nlminb <- tryCatch(
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
        strict_nlminb
      ) &&
      is.finite(
        strict_nlminb$objective
      )
    ) {

      candidates$strict_nlminb <- list(
        par = strict_nlminb$par,
        objective =
          strict_nlminb$objective,
        convergence =
          strict_nlminb$convergence,
        message =
          strict_nlminb$message
      )
    }

    strict_start <- if (
      !is.null(
        strict_nlminb
      ) &&
      is.finite(
        strict_nlminb$objective
      )
    ) {
      strict_nlminb$par
    } else {
      q_hat
    }

    strict_lbfgsb <- tryCatch(
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
        strict_lbfgsb
      ) &&
      is.finite(
        strict_lbfgsb$value
      )
    ) {

      candidates$strict_L_BFGS_B <- list(
        par = strict_lbfgsb$par,
        objective =
          strict_lbfgsb$value,
        convergence =
          strict_lbfgsb$convergence,
        message =
          strict_lbfgsb$message
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
  }

  gradient_metric_at <- function(
      q
  ) {

    g <- tryCatch(
      bounded_gradient(
        fn = nll,
        x = q,
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

    if (
      !all(
        is.finite(g)
      )
    ) {
      return(
        Inf
      )
    }

    z <- abs(
      g
    )

    if (q[1] <= 1e-6) {
      z[1] <- max(
        -g[1],
        0
      )
    }

    if (q[2] <= 1e-6) {
      z[2] <- max(
        -g[2],
        0
      )
    }

    max(
      z
    )
  }

  post_strict_gradient_metric <-
    gradient_metric_at(
      q_hat
    )

  if (
    !is.finite(
      post_strict_gradient_metric
    ) ||
    post_strict_gradient_metric >=
    0.01
  ) {

    rescue_pool <- list(
      current = fit
    )

    ultra_original_start <- q_hat

    ultra_nlminb <- tryCatch(
      nlminb(
        start = q_hat,
        objective = nll,
        lower = lower,
        upper = upper,
        control = list(
          eval.max = 8000,
          iter.max = 5000,
          rel.tol = 1e-13,
          x.tol = 1e-13
        )
      ),
      error = function(e) {
        NULL
      }
    )

    if (
      !is.null(
        ultra_nlminb
      ) &&
      is.finite(
        ultra_nlminb$objective
      )
    ) {

      rescue_pool$ultra_nlminb <- list(
        par = ultra_nlminb$par,
        objective =
          ultra_nlminb$objective,
        convergence =
          ultra_nlminb$convergence,
        message =
          ultra_nlminb$message
      )
    }

    ultra_start <- if (
      !is.null(
        ultra_nlminb
      ) &&
      is.finite(
        ultra_nlminb$objective
      )
    ) {
      ultra_nlminb$par
    } else {
      q_hat
    }

    ultra_bfgs <- tryCatch(
      optim(
        par = ultra_start,
        fn = nll,
        method = "BFGS",
        control = list(
          maxit = 5000,
          reltol = 1e-13
        )
      ),
      error = function(e) {
        NULL
      }
    )

    if (
      !is.null(
        ultra_bfgs
      ) &&
      is.finite(
        ultra_bfgs$value
      ) &&
      all(
        ultra_bfgs$par >=
          lower
      ) &&
      all(
        ultra_bfgs$par <=
          upper
      )
    ) {

      rescue_pool$ultra_BFGS <- list(
        par = ultra_bfgs$par,
        objective =
          ultra_bfgs$value,
        convergence =
          ultra_bfgs$convergence,
        message =
          ultra_bfgs$message
      )
    }

    if (
      !is.null(
        ultra_bfgs
      ) &&
      is.finite(
        ultra_bfgs$value
      ) &&
      all(
        ultra_bfgs$par >=
          lower
      ) &&
      all(
        ultra_bfgs$par <=
          upper
      )
    ) {

      ultra_bfgs_nlminb <- tryCatch(
        nlminb(
          start = ultra_bfgs$par,
          objective = nll,
          lower = lower,
          upper = upper,
          control = list(
            eval.max = 8000,
            iter.max = 5000,
            rel.tol = 1e-13,
            x.tol = 1e-13
          )
        ),
        error = function(e) {
          NULL
        }
      )

      if (
        !is.null(
          ultra_bfgs_nlminb
        ) &&
        is.finite(
          ultra_bfgs_nlminb$objective
        )
      ) {

        rescue_pool$ultra_BFGS_nlminb <- list(
          par =
            ultra_bfgs_nlminb$par,
          objective =
            ultra_bfgs_nlminb$objective,
          convergence =
            ultra_bfgs_nlminb$convergence,
          message =
            ultra_bfgs_nlminb$message
        )
      }
    }

    ultra_direct_bfgs <- tryCatch(
      optim(
        par = ultra_original_start,
        fn = nll,
        method = "BFGS",
        control = list(
          maxit = 5000,
          reltol = 1e-13
        )
      ),
      error = function(e) {
        NULL
      }
    )

    if (
      !is.null(
        ultra_direct_bfgs
      ) &&
      is.finite(
        ultra_direct_bfgs$value
      ) &&
      all(
        ultra_direct_bfgs$par >=
          lower
      ) &&
      all(
        ultra_direct_bfgs$par <=
          upper
      )
    ) {

      rescue_pool$ultra_direct_BFGS <- list(
        par =
          ultra_direct_bfgs$par,
        objective =
          ultra_direct_bfgs$value,
        convergence =
          ultra_direct_bfgs$convergence,
        message =
          ultra_direct_bfgs$message
      )

      ultra_direct_bfgs_nlminb <- tryCatch(
        nlminb(
          start = ultra_direct_bfgs$par,
          objective = nll,
          lower = lower,
          upper = upper,
          control = list(
            eval.max = 8000,
            iter.max = 5000,
            rel.tol = 1e-13,
            x.tol = 1e-13
          )
        ),
        error = function(e) {
          NULL
        }
      )

      if (
        !is.null(
          ultra_direct_bfgs_nlminb
        ) &&
        is.finite(
          ultra_direct_bfgs_nlminb$objective
        )
      ) {

        rescue_pool$ultra_direct_BFGS_nlminb <- list(
          par =
            ultra_direct_bfgs_nlminb$par,
          objective =
            ultra_direct_bfgs_nlminb$objective,
          convergence =
            ultra_direct_bfgs_nlminb$convergence,
          message =
            ultra_direct_bfgs_nlminb$message
        )
      }
    }

    pool_objective <- vapply(
      rescue_pool,
      function(z) {
        z$objective
      },
      numeric(1)
    )

    pool_convergence <- vapply(
      rescue_pool,
      function(z) {
        z$convergence
      },
      integer(1)
    )

    pool_gradient <- vapply(
      rescue_pool,
      function(z) {
        gradient_metric_at(
          z$par
        )
      },
      numeric(1)
    )

    good <- which(
      pool_convergence ==
        0L &
        is.finite(
          pool_objective
        )
    )

    if (length(good) == 0L) {

      good <- which(
        is.finite(
          pool_objective
        )
      )
    }

    best_objective <- min(
      pool_objective[
        good
      ]
    )

    eligible <- good[
      pool_objective[
        good
      ] <=
        best_objective +
        1e-6
    ]

    chosen <- eligible[
      which.min(
        pool_gradient[
          eligible
        ]
      )
    ]

    fit <- rescue_pool[[
      chosen
    ]]

    selected_method <- names(
      rescue_pool
    )[
      chosen
    ]

    q_hat <- fit$par
  }

  estimate <- c(
    theta = q_hat[1],
    phi = q_hat[2],
    alpha = q_hat[3],
    lambda = exp(
      q_hat[4]
    ),
    rho = exp(
      q_hat[5]
    )
  )

  distance_lower <-
    q_hat -
    lower

  distance_upper <-
    upper -
    q_hat

  theta_boundary <-
    q_hat[1] <=
    1e-6

  phi_boundary <-
    q_hat[2] <=
    1e-6

  other_boundary_hit <- any(
    c(
      distance_upper[
        1:2
      ] <
        1e-5,
      distance_lower[
        3:5
      ] <
        1e-5,
      distance_upper[
        3:5
      ] <
        1e-5
    )
  )

  boundary_hit <-
    theta_boundary ||
    phi_boundary ||
    other_boundary_hit

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

  theta_gradient <- if (
    is.finite(
      gradient[1]
    )
  ) {
    gradient[1]
  } else {
    NA_real_
  }

  phi_gradient <- if (
    is.finite(
      gradient[2]
    )
  ) {
    gradient[2]
  } else {
    NA_real_
  }

  natural_boundary_index <- c(
    theta_boundary,
    phi_boundary,
    FALSE,
    FALSE,
    FALSE
  )

  interior_gradient <-
    gradient[
      !natural_boundary_index
    ]

  max_abs_gradient_interior <- if (
    all(
      is.finite(
        interior_gradient
      )
    )
  ) {
    max(
      abs(
        interior_gradient
      )
    )
  } else {
    NA_real_
  }

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

  gradient_metric <- if (
    all(
      is.finite(
        gradient
      )
    )
  ) {

    metric_components <- abs(
      gradient
    )

    if (theta_boundary) {
      metric_components[1] <- max(
        -gradient[1],
        0
      )
    }

    if (phi_boundary) {
      metric_components[2] <- max(
        -gradient[2],
        0
      )
    }

    max(
      metric_components
    )

  } else {

    NA_real_
  }

  gradient_ok <-
    is.finite(
      gradient_metric
    ) &&
    gradient_metric <
    0.01

  asymmetric_details <- nll(
    q_hat,
    diagnostics = TRUE
  )

  nll_reference <- make_nested_nll(
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
    !other_boundary_hit &&
    is.finite(
      fit$objective
    ) &&
    gradient_ok &&
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
      theta_boundary =
        theta_boundary,
      phi_boundary =
        phi_boundary,
      theta_gradient =
        theta_gradient,
      phi_gradient =
        phi_gradient,
      gradient_metric =
        gradient_metric,
      max_abs_gradient =
        max_abs_gradient,
      max_abs_gradient_interior =
        max_abs_gradient_interior,
      gradient_ok =
        gradient_ok,
      pre_refinement_gradient_metric =
        pre_refinement_gradient_metric,
      strict_refinement_attempted =
        strict_refinement_attempted,
      asymmetric_logLik =
        loglik_asymmetric,
      reference_logLik =
        loglik_reference,
      quadrature_difference =
        quadrature_difference,
      quadrature_nodes = 64L,
      quadrature_log_drop = 20,
      v_integral_count =
        asymmetric_details$v_integral_count,
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
      mean_v =
        sim$mean_v,
      var_v =
        sim$var_v,
      mean_w =
        sim$mean_w,
      var_w =
        sim$var_w,
      mean_product_frailty =
        sim$mean_product_frailty,
      var_product_frailty =
        sim$var_product_frailty
    )
  )
}
