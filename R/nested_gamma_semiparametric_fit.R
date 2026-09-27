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
