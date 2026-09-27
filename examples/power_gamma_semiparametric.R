source("R/simulation_utils.R")
source("R/fit_utils.R")
source("R/power_gamma_simulation.R")
source("R/power_gamma_parametric_fit.R")
source("R/power_gamma_semiparametric_fit.R")

set.seed(1234)

sim <- simulate_power_gamma_data(
  n_pairs = 1000,
  theta = 1,
  gamma = 2,
  alpha = -0.9,
  lambda = 0.1,
  rho = 2,
  beta0 = 1,
  beta1 = -0.1,
  entry_max = 3,
  delayed_entry_prob = 0.5,
  censoring_max = 30,
  mz_proportion = 0.5,
  visit_step = 1
)

fit <- fit_power_gamma_semiparametric(
  sim
)

c(
  theta = fit$theta,
  gamma = fit$gamma,
  alpha = fit$alpha
)

fit$fit_ok
fit$logLik
head(
  data.frame(
    event_time = fit$event_times,
    delta = fit$delta
  )
)
