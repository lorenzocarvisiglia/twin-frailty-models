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
