# Predict an Observed Clinical Effect from Two Observed Surrogate Effects

Generates equally weighted predictive draws for a held-out trial's
observed clinical treatment-effect estimate using an already fitted
historical posterior. This function does not fit a historical model and
does not use the held-out clinical estimate to construct predictions.

## Usage

``` r
clinical_predictive_distribution(
  mcmc_dat,
  sample_dat,
  n_draws = NULL,
  seed = NULL,
  diffuse = FALSE,
  diffuse_se = 100,
  intercept0 = FALSE
)
```

## Arguments

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

- n_draws:

  Number of predictive draws to generate. Defaults to the number of rows
  in `mcmc_dat`.

- seed:

  Optional positive integer for reproducible simulation. When supplied,
  the caller's random-number state is restored on exit.

- diffuse:

  Logical; replace the historical surrogate-effect distribution with a
  diffuse prior, as in
  [`MBE()`](https://hyejung0.github.io/MBE/reference/MBE.md). Defaults
  to `FALSE` (full carryover).

- diffuse_se:

  Positive standard deviation of the diffuse surrogate prior.

- intercept0:

  Logical; center the clinical-intercept draws at zero, as in
  [`MBE()`](https://hyejung0.github.io/MBE/reference/MBE.md). Defaults
  to `FALSE`. An absent clinical-intercept column in a fixed-intercept
  fit is treated as zero.

## Value

A numeric vector of `n_draws` equally weighted predictive draws.

## Details

For each historical draw, the joint distribution of the observed
clinical and surrogate estimates is multivariate normal, with the model
covariance plus the within-trial sampling covariance. The function
conditions this joint distribution on the two observed surrogate
estimates using normal conditional-distribution formulas. Historical
parameter draws are reweighted by the marginal likelihood of those
surrogate estimates, then sampled to form the predictive mixture.

This includes residual clinical heterogeneity, clinical sampling
variance, surrogate measurement uncertainty, posterior parameter
uncertainty, and the supplied within-trial sampling correlations. The
clinical estimate itself is never used for weighting or conditioning.
The returned draws target an observed clinical estimate, not a latent
true clinical effect.

This joint conditioning also accounts for clinical-surrogate sampling
correlations and updates historical-draw weights. It differs from the
earlier development assessment's unweighted surrogate-update simulation,
so new predictive intervals need not reproduce archived coverage
results.

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
predictive <- clinical_predictive_distribution(
  historical_posterior[1:100, ], held_out, seed = 2026
)
stats::quantile(predictive, c(0.025, 0.975))
#>       2.5%      97.5% 
#> -0.3798877  0.1683905 
```
