# Check Capture of an Observed Clinical Effect in Predictive Intervals

Calls
[`clinical_predictive_distribution()`](https://hyejung0.github.io/MBE/reference/clinical_predictive_distribution.md)
using an already fitted historical posterior, then checks whether the
observed clinical estimate falls within user-specified predictive
intervals. Historical models are not refitted. The observed clinical
estimate is used only for capture.

## Usage

``` r
tolerance_interval_coverage(
  observed,
  mcmc_dat,
  sample_dat,
  probs = c(0.025, 0.975),
  ...
)
```

## Arguments

- observed:

  One finite numeric observed clinical treatment-effect estimate.

- mcmc_dat:

  Historical posterior draws in the format accepted by
  [`MBE()`](https://hyejung0.github.io/MBE/reference/MBE.md). For
  leave-one-trial-out assessment, fit these draws without the trial to
  be predicted, using
  [`historical_model_fit_2surrogates()`](https://hyejung0.github.io/MBE/reference/historical_model_fit_2surrogates.md).
  Fixed- and random-intercept fits and both uncertainty priors are
  supported.

- sample_dat:

  A named list containing `ClnSE`, `Sur1Est`, `Sur1SE`, `Sur2Est`,
  `Sur2SE`, `R1Clin`, `R2Clin`, and `R12` for one held-out trial. An
  optional `ClnEst` entry is ignored. These are the same sampling inputs
  used by [`MBE()`](https://hyejung0.github.io/MBE/reference/MBE.md),
  except the clinical estimate is not needed.

- probs:

  Lower and upper quantile probabilities: either a numeric vector of
  length two, such as `c(0.025, 0.975)`, or a numeric two-column matrix
  (or data frame) with one row per interval. Column 1 is the lower
  probability and column 2 the upper probability. Each pair must satisfy
  `0 <= lower < upper <= 1`. Asymmetric intervals are allowed.

- ...:

  Additional arguments to
  [`clinical_predictive_distribution()`](https://hyejung0.github.io/MBE/reference/clinical_predictive_distribution.md),
  including `n_draws`, `seed`, `diffuse`, `diffuse_se`, and
  `intercept0`.

## Value

A data frame with one row per interval: `lower_prob`, `upper_prob`,
`level` (100 times their difference), `observed`, interval bounds
`lower` and `upper`, and logical `covered`. For one trial, `covered` is
a capture indicator rather than an aggregate coverage rate.

## Details

Predictive draws are generated once and shared across all intervals in a
call. Bounds use
[`stats::quantile()`](https://rdrr.io/r/stats/quantile.html) with
`type = 7`. Capture is inclusive:
`lower <= observed && observed <= upper`. Missing or infinite inputs are
rejected. Results follow the supplied row order. Specify `seed` for
reproducibility; finite Monte Carlo draws can affect near-boundary
results.

The name follows the manuscript's tolerance-interval terminology. These
are posterior predictive intervals, not classical tolerance intervals
with a separate confidence level. Repeat this function across held-out
trials and average `covered` for each probability pair to obtain
aggregate coverage.

Earlier development code used empirical percentiles and a different
predictive simulation. Archived `loo_assessment_by_trial` indicators
retain that earlier calculation; see
[`clinical_predictive_distribution()`](https://hyejung0.github.io/MBE/reference/clinical_predictive_distribution.md)
for the current conditioning and weighting method.

## Examples

``` r
data("historical_posterior")
data("trial_sim_dat")
held_out <- list(
  ClnSE = trial_sim_dat$CE_se[1],
  Sur1Est = trial_sim_dat$Sur1_est[1],
  Sur1SE = trial_sim_dat$Sur1_se[1],
  Sur2Est = trial_sim_dat$Sur2_est[1],
  Sur2SE = trial_sim_dat$Sur2_se[1],
  R1Clin = trial_sim_dat$Cor_CE_Sur1[1],
  R2Clin = trial_sim_dat$Cor_CE_Sur2[1],
  R12 = trial_sim_dat$Cor_Sur1_Sur2[1]
)
tolerance_interval_coverage(
  observed = trial_sim_dat$CE_est[1],
  mcmc_dat = historical_posterior[1:100, ], sample_dat = held_out,
  probs = rbind(c(0.025, 0.975), c(0.05, 0.95), c(0.075, 0.925)),
  seed = 2026
)
#>   lower_prob upper_prob level   observed      lower      upper covered
#> 1      0.025      0.975    95 -0.1337531 -0.3798877 0.16839047    TRUE
#> 2      0.050      0.950    90 -0.1337531 -0.3472340 0.13095301    TRUE
#> 3      0.075      0.925    85 -0.1337531 -0.3010150 0.07382694    TRUE
```
