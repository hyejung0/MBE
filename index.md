# Multi-component Bayesian Endpoints (MBEs)

This repository contains all relevant R codes that was used in the
Multi-Component Bayesian Endpoints for Randomized Clinical Trials paper
produced by Lee et. al., (2026).

This package implements an MBE with exactly one clinical endpoint and
two surrogate endpoints. The current CKD implementation uses chronic and
acute eGFR slopes as the two surrogates. A Bayesian framework
incorporates historical randomized clinical trials (RCTs) and propagates
uncertainty when estimating the treatment effect on the clinical
endpoint from both surrogate effects.

Detailed explanation of mapping the manuscript to the code can be found
in this document: [Mapping the manuscript to the
code](https://hyejung0.github.io/MBE/articles/MBE-manuscript-to-code.md)

## Installation

You can install the development version of MBE from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pak("hyejung0/MBE")
```

## Example

The package includes compatible historical posterior draws and simulated
data. The input for a new trial always contains the clinical endpoint,
both surrogate endpoints, and the three pairwise correlations:

``` r

library(MBE)

data("historical_posterior")
data("trial_sim_dat")

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

set.seed(1)
fit <- MBE(historical_posterior, new_trial)
fit$post_mean
```
