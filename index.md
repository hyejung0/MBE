# Multi-component Bayesian Endpoints (MBEs)

This repository contains the R package accompanying the *Multi-Component
Bayesian Endpoints for Randomized Clinical Trials* manuscript by Lee et
al. (2026).

This package implements an MBE with exactly one clinical endpoint and
two surrogate endpoints. The current CKD implementation uses chronic and
acute eGFR slopes as the two surrogates. A Bayesian framework
incorporates historical randomized clinical trials (RCTs) and propagates
uncertainty when estimating the treatment effect on the clinical
endpoint from both surrogate effects.

See [Mapping the MBE Manuscript to the
Package](https://hyejung0.github.io/MBE/articles/MBE-manuscript-to-code.html)
for the complete workflow.

## Installation

You can install the development version of MBE from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pak("hyejung0/MBE")
```

Fitting a historical model requires a working CmdStan installation. If
CmdStan has not already been installed, run:

``` r

cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan()
```

## Example

The package includes compatible historical posterior draws and simulated
data. The input for a new trial always contains the clinical endpoint,
both surrogate endpoints, and the three pairwise correlations:

``` r

library(MBE)

data("historical_posterior")
data("trial_sim_dat")

#Construct a list for the new trial's observed data, using the first row of the simulated dataset.
new_trial <- list(
  ClnEst = trial_sim_dat$CE_est[1],
  ClnSE = trial_sim_dat$CE_se[1],
  Sur1Est = trial_sim_dat$Sur1_est[1],
  Sur1SE = trial_sim_dat$Sur1_se[1],
  Sur2Est = trial_sim_dat$Sur2_est[1],
  Sur2SE = trial_sim_dat$Sur2_se[1],
  R1Clin = trial_sim_dat$Cor_CE_Sur1[1],
  R2Clin = trial_sim_dat$Cor_CE_Sur2[1],
  R12 = trial_sim_dat$Cor_Sur1_Sur2[1]
)

#Estimate the posterior distribution of the MBE for the new trial using the historical posterior draws.
set.seed(1)
fit <- MBE(historical_posterior, new_trial)

#Mean of the posterior distribution for the clinical endpoint, surrogate 1, and surrogate 2.
fit$post_mean
#>                  [,1]
#> clinical   -0.1013294
#> surrogate1  0.6837417
#> surrogate2 -5.4760921

#Corresponding variance matrix
fit$post_var
#>               clinical   surrogate1   surrogate2
#> clinical    0.01196417 -0.019592518 -0.091338433
#> surrogate1 -0.01959252  0.060139164 -0.009146439
#> surrogate2 -0.09133843 -0.009146439  3.981621682
```

A detailed example is provided in the package vignette [Simulation
Example](https://hyejung0.github.io/MBE/articles/simulation-example.html).
