source("R/simulation_utils.R")
source("R/fit_utils.R")
source("R/correlated_lognormal_simulation.R")
source("R/correlated_lognormal_fit.R")

set.seed(1234)

sim <- simulate_correlated_lognormal_data(
  n_pairs = 1000,
  sigma2 = 1,
  gamma = 0.5,
  alpha = -0.9,
  lambda = 0.1,
  rho = 2,
  beta0 = 1,
  beta1 = -0.01,
  entry_max = 3,
  delayed_entry_prob = 0.5,
  censoring_max = 30,
  mz_proportion = 0.5,
  visit_step = 1
)

cat("Generated pairs:", sim$n_pairs_generated, "\n")
cat("Observed pairs:", sim$n_pairs_observed, "\n")
cat("Observed DZ pairs:", sim$n_dz_observed, "\n")
cat("Observed MZ pairs:", sim$n_mz_observed, "\n")
cat("Observed subjects:", sim$n_subjects_observed, "\n")
cat("Events:", sim$n_events, "\n")
cat("Event proportion:", sim$event_proportion, "\n")
cat("Mean W:", sim$mean_w, "\n")
cat("Variance W:", sim$var_w, "\n")
cat("Frailty mean:", sim$mean_frailty, "\n")
cat("Frailty variance:", sim$var_frailty, "\n")

head(sim$subjects)
head(sim$long)

fit <- fit_correlated_lognormal(
  sim
)

fit$estimates
fit$fit_ok
fit$diagnostics
