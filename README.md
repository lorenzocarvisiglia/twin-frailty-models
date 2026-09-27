# Twin frailty models

R code for simulation and estimation of frailty survival models for twin data with delayed entry, right censoring, and longitudinal time-varying covariates.

The repository implements three frailty structures:

- power gamma frailty
- nested gamma frailty
- correlated lognormal frailty

## Repository structure

```text
R/
├── simulation_utils.R
├── fit_utils.R
├── power_gamma_simulation.R
├── power_gamma_fit.R
├── nested_gamma_simulation.R
├── nested_gamma_fit.R
├── correlated_lognormal_simulation.R
└── correlated_lognormal_fit.R

examples/
├── power_gamma_simulation.R
├── nested_gamma_simulation.R
└── correlated_lognormal_simulation.R
