normal_quadrature_rule <- function(
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

    off <- sqrt(
      j / 2
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

  nodes <- eig$values[
    ord
  ]

  weights <- eig$vectors[
    1,
    ord
  ]^2

  list(
    z = sqrt(2) *
      nodes,
    log_weights = log(
      pmax(
        weights,
        .Machine$double.xmin
      )
    )
  )
}


row_log_sum_exp <- function(
    x
) {

  m <- apply(
    x,
    1,
    max
  )

  m +
    log(
      rowSums(
        exp(
          sweep(
            x,
            1,
            m,
            "-"
          )
        )
      )
    )
}


prepare_correlated_lognormal_fit <- function(
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


make_correlated_lognormal_nll <- function(
    prep,
    n_nodes = 25L
) {

  rule <- normal_quadrature_rule(
    n_nodes
  )

  z <- rule$z
  log_weights <- rule$log_weights

  z1_2d <- rep(
    z,
    each = n_nodes
  )

  z2_2d <- rep(
    z,
    times = n_nodes
  )

  log_weights_2d <-
    rep(
      log_weights,
      each = n_nodes
    ) +
    rep(
      log_weights,
      times = n_nodes
    )

  function(q) {

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

    sigma2 <- q[1]
    gamma <- q[2]
    alpha <- q[3]

    lambda <- exp(
      q[4]
    )

    rho <- exp(
      q[5]
    )

    if (
      sigma2 < 0 ||
      abs(gamma) > 1
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

    sigma <- sqrt(
      sigma2
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

      d <-
        prep$subject_status[
          s1
        ] +
        prep$subject_status[
          s2
        ]

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

      w <- sigma *
        z

      exp_w <- exp(
        w
      )

      log_num <-
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

      log_num <- sweep(
        log_num,
        2,
        log_weights,
        "+"
      )

      log_den <-
        -outer(
          A,
          exp_w,
          "*"
        )

      log_den <- sweep(
        log_den,
        2,
        log_weights,
        "+"
      )

      loglik_mz <- sum(
        row_log_sum_exp(
          log_num
        ) -
          row_log_sum_exp(
            log_den
          )
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

      w1 <- sigma *
        z1_2d

      if (gamma == -1) {

        w2 <- -w1

      } else if (gamma == 1) {

        w2 <- w1

      } else {

        w2 <- sigma *
          (
            gamma *
              z1_2d +
              sqrt(
                1 -
                  gamma^2
              ) *
              z2_2d
          )
      }

      exp_w1 <- exp(
        w1
      )

      exp_w2 <- exp(
        w2
      )

      log_num <-
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

      log_num <- sweep(
        log_num,
        2,
        log_weights_2d,
        "+"
      )

      log_den <-
        -outer(
          A1,
          exp_w1,
          "*"
        ) -
        outer(
          A2,
          exp_w2,
          "*"
        )

      log_den <- sweep(
        log_den,
        2,
        log_weights_2d,
        "+"
      )

      loglik_dz <- sum(
        row_log_sum_exp(
          log_num
        ) -
          row_log_sum_exp(
            log_den
          )
      )
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

    -loglik
  }
}


fit_correlated_lognormal <- function(
    sim
) {

  prep <- prepare_correlated_lognormal_fit(
    sim
  )

  nll <- make_correlated_lognormal_nll(
    prep,
    n_nodes = 25L
  )

  lower <- c(
    0,
    -1,
    -10,
    log(1e-5),
    log(0.3)
  )

  upper <- c(
    10,
    1,
    10,
    log(2),
    log(8)
  )

  start1 <- c(
    0.8,
    0.3,
    0,
    log(0.02),
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
      1.5,
      0,
      -1,
      log(0.01),
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
      "all correlated lognormal optimizers failed"
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

  refined <- tryCatch(
    nlminb(
      start = fit$par,
      objective = nll,
      lower = lower,
      upper = upper,
      control = list(
        eval.max = 5000,
        iter.max = 3000,
        rel.tol = 1e-12,
        x.tol = 1e-12
      )
    ),
    error = function(e) {
      NULL
    }
  )

  if (
    !is.null(refined) &&
    refined$convergence ==
      0L &&
    is.finite(
      refined$objective
    ) &&
    refined$objective <=
      fit$objective +
      1e-4
  ) {

    fit <- list(
      par = refined$par,
      objective = refined$objective,
      convergence = refined$convergence,
      message = refined$message
    )

    selected_method <- paste0(
      selected_method,
      "+strict_nlminb"
    )
  }

  q_pre_high <- fit$par

  nll_50 <- make_correlated_lognormal_nll(
    prep,
    n_nodes = 50L
  )

  nll_60 <- make_correlated_lognormal_nll(
    prep,
    n_nodes = 60L
  )

  high_accuracy_gradient_status <- function(
      q,
      fn
  ) {

    distance_lower_local <-
      q -
      lower

    distance_upper_local <-
      upper -
      q

    sigma2_boundary_local <-
      q[1] <=
      1e-6

    gamma_lower_boundary_local <-
      distance_lower_local[2] <
      1e-5

    gamma_upper_boundary_local <-
      distance_upper_local[2] <
      1e-5

    gradient_local <- tryCatch(
      bounded_gradient(
        fn = fn,
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
        is.finite(
          gradient_local
        )
      )
    ) {
      return(
        list(
          metric = Inf,
          gradient = gradient_local
        )
      )
    }

    metric_components <- numeric()

    if (sigma2_boundary_local) {

      metric_components <- c(
        metric_components,
        max(
          -gradient_local[1],
          0
        )
      )

    } else {

      metric_components <- c(
        metric_components,
        abs(
          gradient_local[1]
        )
      )
    }

    if (!sigma2_boundary_local) {

      if (gamma_lower_boundary_local) {

        metric_components <- c(
          metric_components,
          max(
            -gradient_local[2],
            0
          )
        )

      } else if (gamma_upper_boundary_local) {

        metric_components <- c(
          metric_components,
          max(
            gradient_local[2],
            0
          )
        )

      } else {

        metric_components <- c(
          metric_components,
          abs(
            gradient_local[2]
          )
        )
      }
    }

    metric_components <- c(
      metric_components,
      abs(
        gradient_local[
          3:5
        ]
      )
    )

    list(
      metric = max(
        metric_components
      ),
      gradient = gradient_local
    )
  }

  high_pre_status <-
    high_accuracy_gradient_status(
      q = q_pre_high,
      fn = nll_60
    )

  high_accuracy_pre_gradient_metric <-
    high_pre_status$metric

  high_accuracy_pre_quad_diff <-
    (
      -nll_50(
        q_pre_high
      )
    ) -
    (
      -nll_60(
        q_pre_high
      )
    )

  high_accuracy_refinement <-
    !is.finite(
      high_accuracy_pre_gradient_metric
    ) ||
    high_accuracy_pre_gradient_metric >=
      0.01 ||
    !is.finite(
      high_accuracy_pre_quad_diff
    ) ||
    abs(
      high_accuracy_pre_quad_diff
    ) >=
      0.001

  if (high_accuracy_refinement) {

    gr_60 <- function(q) {
      bounded_gradient(
        fn = nll_60,
        x = q,
        lower = lower,
        upper = upper,
        rel_step = 1e-4
      )
    }

    high_candidates <- list(
      start25 = list(
        par = q_pre_high,
        objective = nll_60(
          q_pre_high
        ),
        convergence = fit$convergence,
        message = "25-node estimate retained"
      )
    )

    high_fit1 <- tryCatch(
      optim(
        par = q_pre_high,
        fn = nll_60,
        gr = gr_60,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(
          maxit = 4000,
          factr = 1e4,
          pgtol = 1e-8
        )
      ),
      error = function(e) {
        NULL
      }
    )

    if (
      !is.null(
        high_fit1
      ) &&
      is.finite(
        high_fit1$value
      )
    ) {

      high_candidates$L_BFGS_B <- list(
        par = high_fit1$par,
        objective = high_fit1$value,
        convergence = high_fit1$convergence,
        message = high_fit1$message
      )
    }

    high_start2 <- if (
      !is.null(
        high_fit1
      ) &&
      is.finite(
        high_fit1$value
      )
    ) {
      high_fit1$par
    } else {
      q_pre_high
    }

    high_fit2 <- tryCatch(
      nlminb(
        start = high_start2,
        objective = nll_60,
        gradient = gr_60,
        lower = lower,
        upper = upper,
        control = list(
          eval.max = 6000,
          iter.max = 4000,
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
        high_fit2
      ) &&
      is.finite(
        high_fit2$objective
      )
    ) {

      high_candidates$nlminb <- list(
        par = high_fit2$par,
        objective = high_fit2$objective,
        convergence = high_fit2$convergence,
        message = high_fit2$message
      )
    }

    if (
      !is.null(
        high_fit2
      ) &&
      is.finite(
        high_fit2$objective
      ) &&
      high_fit2$convergence !=
        0L
    ) {

      high_fit3 <- tryCatch(
        optim(
          par = high_fit2$par,
          fn = nll_60,
          gr = gr_60,
          method = "L-BFGS-B",
          lower = lower,
          upper = upper,
          control = list(
            maxit = 5000,
            factr = 1e4,
            pgtol = 1e-8
          )
        ),
        error = function(e) {
          NULL
        }
      )

      if (
        !is.null(
          high_fit3
        ) &&
        is.finite(
          high_fit3$value
        )
      ) {

        high_candidates$rescue_L_BFGS_B <- list(
          par = high_fit3$par,
          objective = high_fit3$value,
          convergence = high_fit3$convergence,
          message = high_fit3$message
        )
      }
    }

    high_objectives <- vapply(
      high_candidates,
      function(x) {
        x$objective
      },
      numeric(1)
    )

    high_metrics <- vapply(
      high_candidates,
      function(x) {
        high_accuracy_gradient_status(
          q = x$par,
          fn = nll_60
        )$metric
      },
      numeric(1)
    )

    high_convergence_zero <- vapply(
      high_candidates,
      function(x) {
        x$convergence ==
          0L
      },
      logical(1)
    )

    high_best_objective <- min(
      high_objectives
    )

    high_eligible <- which(
      high_objectives <=
        high_best_objective +
        1e-6
    )

    high_eligible_converged <-
      high_eligible[
        high_convergence_zero[
          high_eligible
        ]
      ]

    if (
      length(
        high_eligible_converged
      ) >
        0L
    ) {
      high_eligible <-
        high_eligible_converged
    }

    high_selected_idx <-
      high_eligible[
        which.min(
          high_metrics[
            high_eligible
          ]
        )
      ]

    high_selected <-
      high_candidates[[
        high_selected_idx
      ]]

    fit <- list(
      par = high_selected$par,
      objective = high_selected$objective,
      convergence = high_selected$convergence,
      message = high_selected$message
    )

    selected_method <- paste0(
      selected_method,
      "+high60_",
      names(
        high_candidates
      )[
        high_selected_idx
      ]
    )
  }

  q_hat <- fit$par

  fit$objective <- nll_60(
    q_hat
  )

  estimate <- c(
    sigma2 = q_hat[1],
    gamma = q_hat[2],
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

  sigma2_boundary <-
    q_hat[1] <=
    1e-6

  gamma_lower_boundary <-
    distance_lower[2] <
    1e-5

  gamma_upper_boundary <-
    distance_upper[2] <
    1e-5

  gamma_boundary <-
    gamma_lower_boundary ||
    gamma_upper_boundary

  other_boundary_hit <- any(
    distance_lower[
      3:5
    ] <
      1e-5 |
      distance_upper[
        3:5
      ] <
      1e-5
  )

  boundary_hit <-
    sigma2_boundary ||
    gamma_boundary ||
    other_boundary_hit

  gradient <- tryCatch(
    bounded_gradient(
      fn = nll_60,
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

  sigma2_gradient <- gradient[1]
  gamma_gradient <- gradient[2]

  interior_indices <- 1:5

  if (sigma2_boundary) {
    interior_indices <- setdiff(
      interior_indices,
      1L
    )
  }

  if (
    gamma_boundary ||
    sigma2_boundary
  ) {
    interior_indices <- setdiff(
      interior_indices,
      2L
    )
  }

  max_abs_gradient_interior <- if (
    length(interior_indices) ==
      0L
  ) {

    0

  } else if (
    all(
      is.finite(
        gradient[
          interior_indices
        ]
      )
    )
  ) {

    max(
      abs(
        gradient[
          interior_indices
        ]
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

  sigma2_kkt_ok <-
    if (sigma2_boundary) {

      is.finite(
        sigma2_gradient
      ) &&
        sigma2_gradient >=
        -0.01

    } else {

      TRUE
    }

  gamma_kkt_ok <-
    if (sigma2_boundary) {

      TRUE

    } else if (gamma_lower_boundary) {

      is.finite(
        gamma_gradient
      ) &&
        gamma_gradient >=
        -0.01

    } else if (gamma_upper_boundary) {

      is.finite(
        gamma_gradient
      ) &&
        gamma_gradient <=
        0.01

    } else {

      is.finite(
        gamma_gradient
      ) &&
        abs(
          gamma_gradient
        ) <
        0.01
    }

  gradient_ok <-
    sigma2_kkt_ok &&
    gamma_kkt_ok &&
    is.finite(
      max_abs_gradient_interior
    ) &&
    max_abs_gradient_interior <
      0.01

  nll_20 <- make_correlated_lognormal_nll(
    prep,
    n_nodes = 20L
  )

  nll_30 <- make_correlated_lognormal_nll(
    prep,
    n_nodes = 30L
  )

  loglik_20 <- -nll_20(
    q_hat
  )

  loglik_30 <- -nll_30(
    q_hat
  )

  quadrature_difference <-
    loglik_20 -
    loglik_30

  loglik_50 <- -nll_50(
    q_hat
  )

  loglik_60 <- -nll_60(
    q_hat
  )

  quadrature_difference_50_60 <-
    loglik_50 -
    loglik_60

  fit_ok <-
    fit$convergence ==
    0L &&
    !other_boundary_hit &&
    is.finite(
      fit$objective
    ) &&
    gradient_ok &&
    is.finite(
      quadrature_difference_50_60
    ) &&
    abs(
      quadrature_difference_50_60
    ) <
    0.001

  list(
    estimates = estimate,
    convergence =
      fit$convergence,
    fit_ok =
      fit_ok,
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
    logLik =
      -fit$objective,
    n_events =
      sim$n_events,
    event_proportion =
      sim$event_proportion,
    diagnostics = list(
      selected_method =
        selected_method,
      high_accuracy_refinement =
        high_accuracy_refinement,
      high_accuracy_pre_gradient_metric =
        high_accuracy_pre_gradient_metric,
      high_accuracy_pre_quad_diff =
        high_accuracy_pre_quad_diff,
      objective =
        fit$objective,
      boundary_hit =
        boundary_hit,
      sigma2_boundary =
        sigma2_boundary,
      gamma_boundary =
        gamma_boundary,
      gamma_identifiable =
        !sigma2_boundary,
      sigma2_gradient =
        sigma2_gradient,
      gamma_gradient =
        gamma_gradient,
      gamma_lower_boundary =
        gamma_lower_boundary,
      gamma_upper_boundary =
        gamma_upper_boundary,
      gamma_kkt_ok =
        gamma_kkt_ok,
      max_abs_gradient =
        max_abs_gradient,
      max_abs_gradient_interior =
        max_abs_gradient_interior,
      gradient_ok =
        gradient_ok,
      quadrature_logLik_20 =
        loglik_20,
      quadrature_logLik_30 =
        loglik_30,
      quadrature_difference =
        quadrature_difference,
      quadrature_logLik_50 =
        loglik_50,
      quadrature_logLik_60 =
        loglik_60,
      quadrature_difference_50_60 =
        quadrature_difference_50_60,
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
      mean_w =
        sim$mean_w,
      var_w =
        sim$var_w,
      mean_frailty =
        sim$mean_frailty,
      var_frailty =
        sim$var_frailty
    )
  )
}
