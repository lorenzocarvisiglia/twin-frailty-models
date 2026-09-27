# Twin frailty models

R code for simulation and estimation of frailty survival models for twin data with delayed entry, right censoring, and longitudinal time-varying covariates.

The repository implements three frailty structures:

- power gamma frailty
- nested gamma frailty
- correlated lognormal frailty

For each frailty structure, the repository provides:

- simulation under a Weibull baseline hazard
- parametric maximum-likelihood estimation with a Weibull baseline
- semiparametric estimation with an unspecified baseline cumulative hazard represented by jumps at the observed event times

Delayed entry is handled by conditioning on survival to the individual entry time. The covariate contribution before entry uses the covariate value at entry, while the observed post-entry time-varying covariate is represented in start-stop form.

## Repository structure

```text
R/
├── simulation_utils.R
├── fit_utils.R
├── power_gamma_simulation.R
├── power_gamma_parametric_fit.R
├── power_gamma_semiparametric_fit.R
├── nested_gamma_simulation.R
├── nested_gamma_parametric_fit.R
├── nested_gamma_semiparametric_fit.R
├── correlated_lognormal_simulation.R
├── correlated_lognormal_parametric_fit.R
└── correlated_lognormal_semiparametric_fit.R

examples/
├── power_gamma_simulation.R
├── power_gamma_semiparametric.R
├── nested_gamma_simulation.R
├── nested_gamma_semiparametric.R
├── correlated_lognormal_simulation.R
└── correlated_lognormal_semiparametric.R
```

## Parametric examples

The existing `*_simulation.R` examples simulate one data set and fit the corresponding Weibull-baseline model.

For example:

```r
source("R/simulation_utils.R")
source("R/fit_utils.R")
source("R/power_gamma_simulation.R")
source("R/power_gamma_parametric_fit.R")

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

fit <- fit_power_gamma(sim)
fit$estimates
```

The corresponding parametric fitting functions are:

```r
fit_power_gamma(sim)
fit_nested_gamma(sim)
fit_correlated_lognormal(sim)
```

## Semiparametric examples

The semiparametric fitters use the corresponding parametric fit for initialization, so both fitting files should be sourced.

The main semiparametric fitting functions are:

```r
fit_power_gamma_semiparametric_profile(sim)
fit_nested_gamma_semiparametric_profile(sim)
fit_correlated_lognormal_semiparametric_profile(sim)
```

Complete examples are provided in the `examples/` directory.

The fitted baseline is returned through:

- `event_times`: observed event times
- `delta`: estimated jumps of the baseline cumulative hazard

The structural frailty parameters and the time-varying covariate coefficient are returned by the model-specific fit object together with convergence, gradient, boundary, and quadrature diagnostics.

## Data representation

Simulation functions return subject-level and start-stop longitudinal survival data. The fitting functions use complete twin pairs and distinguish MZ and DZ contributions according to the selected frailty structure.

The common utilities are contained in:

- `R/simulation_utils.R` for data generation and start-stop construction
- `R/fit_utils.R` for preparation of fitting data and numerical utilities
