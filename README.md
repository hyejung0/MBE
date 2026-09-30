
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Multi-component Bayesian Endpoints (MBEs)

<!-- badges: start -->

[![R-CMD-check](https://github.com/hyejung0/MBE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hyejung0/MBE/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

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
#>                  [,1]
#> clinical   -0.1013294
#> surrogate1  0.6837417
#> surrogate2 -5.4760921
cat(sprintf(
  "%s: %.3f\n",
  names(fit$importance_diagnostics),
  fit$importance_diagnostics
))
#> effective_sample_size: 3943.574
#>  relative_effective_sample_size: 0.986
#>  maximum_normalized_weight: 0.000
```

## Historical-model assessment

The package provides a restartable workflow for fitting one historical
model per held-out trial and calculating LOO-CV RMSE and 90% and 95%
predictive interval coverage. See `?fit_loo_historical_models` and
`?loo_cv_model_assessment` for the computational workflow.

The completed 66-trial assessment is included as two compact datasets:

``` r
data("loo_assessment_summary", package = "MBE")
data("loo_assessment_by_trial", package = "MBE")

loo_assessment_summary
#>   n_trials loo_cv_rmse coverage_95 coverage_90
#> 1       66   0.2385382    0.969697   0.9090909
head(loo_assessment_by_trial)
#>   trial_id observed_clinical posterior_mean_clinical squared_error
#> 1        1        -0.1337531             -0.09656727   0.001382784
#> 2        2         0.3126354              0.06360662   0.062015334
#> 3        3        -0.7966332             -0.62698427   0.028780758
#> 4        4        -0.3433848             -0.44253095   0.009829956
#> 5        5        -0.1624861             -0.06922007   0.008698544
#> 6        6        -0.5146880             -0.44928852   0.004277089
#>   predictive_percentile predictive_lower_95 predictive_lower_90
#> 1              0.366625          -0.4208025          -0.3664697
#> 2              0.956875          -0.5601490          -0.4838673
#> 3              0.121375          -0.9534428          -0.8945087
#> 4              0.579625          -1.3476991          -1.2068675
#> 5              0.406125          -0.9828023          -0.8245005
#> 6              0.254750          -0.7040234          -0.6563685
#>   predictive_upper_90 predictive_upper_95 covered_95 covered_90
#> 1           0.2232711           0.2789574       TRUE       TRUE
#> 2           0.2948599           0.3647947       TRUE      FALSE
#> 3          -0.2279492          -0.1651125       TRUE       TRUE
#> 4           0.3414936           0.4838808       TRUE       TRUE
#> 5           0.6862105           0.8377323       TRUE       TRUE
#> 6          -0.1831045          -0.1363449       TRUE       TRUE
```
