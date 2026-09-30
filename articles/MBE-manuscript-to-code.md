# Mapping the MBE Manuscript to the Package

## Overview

The Multi-Component Bayesian Endpoint (MBE) workflow described by Lee et
al. (2026) has two stages:

1.  Fit a trial-level Bayesian meta-analytic model to historical
    randomized clinical trials.
2.  Use draws from that historical model to update the treatment effects
    for a new trial.

The package implementation is deliberately restricted to one clinical
endpoint and exactly two surrogate endpoints. In the chronic kidney
disease example, the two surrogates are treatment effects on chronic and
acute eGFR slopes.

## Input data

[`historical_model_fit_2surrogates()`](https://hyejung0.github.io/MBE/reference/historical_model_fit_2surrogates.md)
expects one row per historical trial and the following nine columns:

| Column | Meaning |
|:---|:---|
| `CE_est`, `CE_se` | Clinical treatment-effect estimate and standard error |
| `Sur1_est`, `Sur1_se` | Surrogate 1 estimate and standard error |
| `Sur2_est`, `Sur2_se` | Surrogate 2 estimate and standard error |
| `Cor_CE_Sur1` | Sampling correlation between the clinical endpoint and surrogate 1 |
| `Cor_CE_Sur2` | Sampling correlation between the clinical endpoint and surrogate 2 |
| `Cor_Sur1_Sur2` | Sampling correlation between the two surrogates |

`trial_sim_dat` provides a shareable simulated example with 66 trials.
The restricted CKD-EPI CT trial data used in the motivating analysis are
not distributed with the package.

## Stage 1: fit the historical model

The historical model can use a fixed or estimated clinical intercept and
either half-normal or inverse-gamma priors for its uncertainty
parameters. The manuscript’s initial specification uses an estimated
(random) intercept and inverse-gamma priors. The first simulated trial
can be held out as follows:

``` r

historical_fit <- historical_model_fit_2surrogates(
  data = trial_sim_dat[-1, ],
  random_intercept = TRUE,
  prior_for_uncertainty = "inverse_gamma",
  nchains = 4,
  ncores = 4,
  niter = 2000,
  nwarmup = 1000,
  seed = 2026,
  adapt_delta = 0.95,
  max_treedepth = 15
)

historical_fit$summary
historical_fit$loo
historical_fit$waic
```

The returned object also contains the CmdStan fit, posterior draws, and
sampler diagnostics. Historical fits should be accepted only after
checking convergence and sampler warnings.

## Stage 2: update the MBE for a new trial

The package includes `historical_posterior`, a compact set of compatible
draws from the random-intercept, inverse-gamma model with trial 1 held
out. The held-out trial is converted to the naming convention used by
[`MBE()`](https://hyejung0.github.io/MBE/reference/MBE.md):

``` r

data("historical_posterior", package = "MBE")

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
mbe_fit <- MBE(
  mcmc_dat = historical_posterior,
  sample_dat = new_trial,
  diffuse = TRUE,
  diffuse_se = 100,
  intercept0 = TRUE
)

mbe_fit$post_mean
#>                  [,1]
#> clinical   -0.1013294
#> surrogate1  0.6837417
#> surrogate2 -5.4760921
mbe_fit$post_var
#>               clinical   surrogate1   surrogate2
#> clinical    0.01196417 -0.019592518 -0.091338433
#> surrogate1 -0.01959252  0.060139164 -0.009146439
#> surrogate2 -0.09133843 -0.009146439  3.981621682
mbe_fit$post_quantiles_clinical
#> quantile_0.025  quantile_0.05   quantile_0.1  quantile_0.25   quantile_0.5 
#>    -0.31485836    -0.28123031    -0.24074706    -0.17402273    -0.09799387 
#>  quantile_0.75   quantile_0.9  quantile_0.95 quantile_0.975 
#>    -0.02681263     0.03707755     0.07905239     0.12341721
mbe_fit$importance_diagnostics
#>          effective_sample_size relative_effective_sample_size 
#>                   3.943574e+03                   9.858936e-01 
#>      maximum_normalized_weight 
#>                   3.262557e-04
```

`post_mean` and `post_var` are the mean and covariance of the updated
posterior mixture. `post_psi0` contains resampled posterior draws in the
fixed order clinical endpoint, surrogate 1, and surrogate 2.
`weight_data$norm_w` contains the normalized importance weights.
`importance_diagnostics` helps identify whether those weights are
concentrated in only a small number of historical draws.

## Historical-model fit assessment

The historical model is assessed by holding out each of the 66 trials in
turn.
[`fit_loo_historical_models()`](https://hyejung0.github.io/MBE/reference/fit_loo_historical_models.md)
fits the random-intercept, inverse-gamma model to the other 65 trials
and saves the posterior columns required for assessment. The batch is
restartable: existing trial files are reused unless `overwrite = TRUE`.

``` r

loo_paths <- fit_loo_historical_models(
  data = trial_sim_dat,
  output_dir = "data-raw/loo-random-inverse-gamma",
  seed = 2026,
  nchains = 4,
  ncores = 4,
  niter = 2000,
  nwarmup = 1000,
  adapt_delta = 0.95,
  max_treedepth = 15
)
```

[`loo_cv_model_assessment()`](https://hyejung0.github.io/MBE/reference/loo_cv_model_assessment.md)
calculates two quantities:

- LOO-CV RMSE compares each held-out clinical estimate with its
  posterior mean from the full-carryover MBE update.
- Tolerance coverage predicts the held-out observed clinical estimate
  after conditioning on only its two observed surrogate estimates.

For posterior draw $`b`$, the predictive clinical estimate has
distribution

``` math
N\left(
  \beta_0^{(b)} + \beta_1^{(b)}\gamma_{01}^{(b)} +
  \beta_2^{(b)}\gamma_{02}^{(b)},
  \lambda_{\theta}^{2(b)} + \sigma^2_{\hat{\theta}_0}
\right).
```

The variance therefore includes both the historical model’s residual
clinical heterogeneity and the held-out trial’s clinical sampling
variance.

``` r

loo_assessment <- loo_cv_model_assessment(
  loo_posteriors = "data-raw/loo-random-inverse-gamma",
  data = trial_sim_dat,
  seed = 2026
)
```

The completed, converged 66-trial assessment is included as two compact
package datasets so users do not need to rerun all 66 Stan models:

``` r

data("loo_assessment_summary", package = "MBE")
data("loo_assessment_by_trial", package = "MBE")

knitr::kable(loo_assessment_summary, digits = 4)
```

| n_trials | loo_cv_rmse | coverage_95 | coverage_90 |
|---------:|------------:|------------:|------------:|
|       66 |      0.2385 |      0.9697 |      0.9091 |

``` r


coverage_counts <- data.frame(
  interval = c("95%", "90%"),
  covered = c(
    sum(loo_assessment_by_trial$covered_95),
    sum(loo_assessment_by_trial$covered_90)
  ),
  total = nrow(loo_assessment_by_trial)
)
knitr::kable(coverage_counts)
```

| interval | covered | total |
|:---------|--------:|------:|
| 95%      |      64 |    66 |
| 90%      |      60 |    66 |

Only these small derived tables are installed with MBE. The complete
leave-one-out posterior files remain in `data-raw/` or suitable external
archival storage.
