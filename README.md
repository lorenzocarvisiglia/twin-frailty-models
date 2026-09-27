# Twin frailty models

R code for frailty survival models with delayed entry and longitudinal
time-varying covariates in twin data.

The repository contains simulation and estimation code for three frailty
structures:

- power gamma frailty model
- nested gamma frailty model
- correlated lognormal frailty model

The methods allow for delayed entry, right censoring, and time-varying
covariates represented in start-stop format.

## Repository structure

- `R/`: reusable functions for simulation, likelihood evaluation, and model fitting
- `scripts/`: reproducible scripts for the simulation experiments
- `examples/`: simple examples illustrating how to use the models

## Simulation studies

The repository includes code for:

- correct-specification experiments
- sample-size experiments
- baseline-hazard misspecification
- semiparametric baseline estimation
- frailty-structure misspecification

## Requirements

The code is written in R and is intended to run on a standard personal
computer. Cluster-specific scripts, parallel execution, and machine-specific
file paths are not required.

## Status

The repository is under active development.
