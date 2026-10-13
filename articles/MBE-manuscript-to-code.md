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

The historical dataset should be either a `data.table` or `data.frame`
with the columns names exactly as described below.:

| Column | Meaning |
|:---|:---|
| `CE_est` | Estimate of the treatment effect on the clinical endpoint |
| `CE_se` | Standard error of the estimated treatment effect on the clinical endpoint |
| `Sur1_est` | Estimate of the treatment effect on the first surrogate endpoint |
| `Sur1_se` | Standard error of the estimated treatment effect on the first surrogate endpoint |
| `Sur2_est` | Estimate of the treatment effect on the second surrogate endpoint |
| `Sur2_se` | Standard error of the estimated treatment effect on the second surrogate endpoint |
| `Cor_CE_Sur1` | Sampling correlation between `CE_est` and `Sur1_est` |
| `Cor_CE_Sur2` | Sampling correlation between `CE_est` and `Sur2_est` |
| `Cor_Sur1_Sur2` | Sampling correlation between `Sur1_est` and `Sur2_est` |

`trial_sim_dat` provides a shareable simulated example with 66 trials.
The restricted CKD-EPI CT trial data used in the motivating analysis are
not distributed with the package.

## Stage 1: fit the historical model

The historical model can use a fixed or estimated clinical intercept and
either half-normal or inverse-gamma priors for its uncertainty
parameters. The manuscript’s initial specification uses an random
intercept \beta_0 and inverse-gamma priors on the variance parameters.
The first simulated trial can be held out as follows:

``` r

historical_fields <- c(
  "CE_est", "CE_se", "Sur1_est", "Sur1_se", "Sur2_est", "Sur2_se",
  "Cor_CE_Sur1", "Cor_CE_Sur2", "Cor_Sur1_Sur2"
)
historical_fit <- historical_model_fit_2surrogates(
  data = trial_sim_dat[-1, ..historical_fields],
  random_intercept = TRUE, # allow for a random intercept.
  prior_for_uncertainty = "inverse_gamma", # prior for the variance parameters.
  nchains = 4, # Number of chains
  ncores = 4, # Number of cores for parallel computation.
  # One core per chain is recommended.
  # If set to 1, then all chains will be run sequentially on a single core.
  niter = 2000, # Number of post-warmup samples to save for each chain
  nwarmup = 1000, # Number of warmup samples to discard for each chain
  seed = 2026,
  adapt_delta = 0.95,
  max_treedepth = 15
)
```

The `historical_fit` contains the CmdStan fit, posterior draws, and
sampler diagnostics. Historical fits should be accepted only after
checking convergence and sampler warnings.

The historical posterior distribution samples can be extracted as
follows:

``` r

# Extract historical posterior draws.
historical_posterior <- data.table(historical_fit$posterior_samples)[, c(
  "alphaCEonSur1Sur2",
  "b1CEonSur1Sur2",
  "b2CEonSur1Sur2",
  "SigSqCEonSur1Sur2",
  "alphaSur1onSur2",
  "bSur1onSur2",
  "SigSqSur1onSur2",
  "muSur2",
  "sigSqSur2"
)]
```

``` r

head(historical_posterior)
#>    alphaCEonSur1Sur2 b1CEonSur1Sur2 b2CEonSur1Sur2 SigSqCEonSur1Sur2
#>                <num>          <num>          <num>             <num>
#> 1:       -0.08905504     -0.3454726    -0.03400061      0.0008365727
#> 2:       -0.09577766     -0.2965096    -0.03875367      0.0004283551
#> 3:       -0.09621942     -0.3070684    -0.03997607      0.0018120535
#> 4:       -0.11465709     -0.2251565    -0.03826521      0.0003920071
#> 5:       -0.11554932     -0.2750167    -0.04394211      0.0001956232
#> 6:       -0.06467173     -0.3250740    -0.03458536      0.0016928317
#>    alphaSur1onSur2 bSur1onSur2 SigSqSur1onSur2    muSur2 sigSqSur2
#>              <num>       <num>           <num>     <num>     <num>
#> 1:       0.5896079  0.02679333       0.2645544 -4.147931  39.55489
#> 2:       0.5630964  0.03016332       0.3067365 -3.264356  33.63631
#> 3:       0.5737157  0.03683171       0.3121466 -3.224047  35.62000
#> 4:       0.6327299  0.04116275       0.1920468 -3.806462  43.19725
#> 5:       0.6613998  0.05244425       0.2292295 -3.374513  44.48720
#> 6:       0.7860056  0.02089685       0.2132951 -2.223918  48.81642
```

The column names of this posterior datasets are different from what’s
used in the manuscript. Here’s the mapping:

| Manuscript name      | Programming name    |
|:---------------------|:--------------------|
| \beta_0              | `alphaCEonSur1Sur2` |
| \beta_1              | `b1CEonSur1Sur2`    |
| \beta_2              | `b2CEonSur1Sur2`    |
| \lambda^2\_\theta    | `SigSqCEonSur1Sur2` |
| \alpha_0             | `alphaSur1onSur2`   |
| \alpha_1             | `bSur1onSur2`       |
| \lambda^2\_\gamma    | `SigSqSur1onSur2`   |
| \mu\_{\gamma 2}      | `muSur2`            |
| \sigma^2\_{\gamma 2} | `sigSqSur2`         |

This historical posterior distribution is used as prior for estimating
the treatment effect on the clinical endpoint for the held-out trial.

## Stage 2: update the MBE for a new trial

The held-out trial (first row of `trial_sim_dat`) is converted to the
naming convention used by
[`MBE()`](https://hyejung0.github.io/MBE/reference/MBE.md), and the
posterior distribution is updated using the
[`MBE()`](https://hyejung0.github.io/MBE/reference/MBE.md) function. The
`diffuse` and `diffuse_se` arguments allow the user to specify a diffuse
prior for the new trial’s surrogate endpoints. The `intercept0` argument
allows the user to shift the mean of the intercept to 0.

``` r

# Convert to the naming convention by MBE() and save as a list.
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
  diffuse = TRUE, # Diffuse prior on surrogate endpoints
  diffuse_se = 100,
  intercept0 = TRUE # Center the intercept to 0.
)

mbe_fit$post_mean
#>                   [,1]
#> clinical   -0.05425464
#> surrogate1  0.67469003
#> surrogate2 -5.04103222
mbe_fit$post_var
#>               clinical   surrogate1   surrogate2
#> clinical    0.01212002 -0.018442933 -0.109817889
#> surrogate1 -0.01844293  0.060857682 -0.005088661
#> surrogate2 -0.10981789 -0.005088661  3.798687086
mbe_fit$post_quantiles_clinical
#> quantile_0.025  quantile_0.05   quantile_0.1  quantile_0.25   quantile_0.5 
#>    -0.27319424    -0.23671208    -0.19579364    -0.12941497    -0.05465971 
#>  quantile_0.75   quantile_0.9  quantile_0.95 quantile_0.975 
#>     0.01970057     0.08628304     0.12682206     0.15939645
mbe_fit$importance_diagnostics
#>          effective_sample_size relative_effective_sample_size 
#>                   7.661266e+03                   9.576583e-01 
#>      maximum_normalized_weight 
#>                   2.141063e-04
```

`post_mean` and `post_var` are the mean and covariance of the updated
posterior mixture. These quantities are calculated using weighted
sampling method shown in Appendix B of the manuscript. `post_psi0`
contains resampled posterior draws in the fixed order clinical endpoint,
surrogate 1, and surrogate 2. `weight_data$norm_w` is the set of
normalized weights, \\ \tilde{w}^{(b)} \\. `importance_diagnostics`
helps identify whether those weights are concentrated in only a small
number of historical draws.
