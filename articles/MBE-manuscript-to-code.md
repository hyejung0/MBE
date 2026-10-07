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
intercept $`\beta_0`$ and inverse-gamma priors on the variance
parameters. The first simulated trial can be held out as follows:

``` r

historical_fields <- c(
  "CE_est", "CE_se", "Sur1_est", "Sur1_se", "Sur2_est", "Sur2_se",
  "Cor_CE_Sur1", "Cor_CE_Sur2", "Cor_Sur1_Sur2"
)
historical_fit <- historical_model_fit_2surrogates(
  data = as.data.frame(trial_sim_dat)[-1, historical_fields],
  random_intercept = TRUE, #allow for a random intercept. Set to FALSE to fix the intercept at zero.
  prior_for_uncertainty = "inverse_gamma", #prior for the variance parameters. If set as "half_normal", then half normal distribution will be modeled on the standard deviation parameters.
  nchains = 4, #Number of chains
  ncores = 4, #Number of cores for parallel computation. One core per chain is recommended. If set to 1, then all chains will be run sequentially on a single core.
  niter = 2000, #Number of post-warmup samples to save for each chain
  nwarmup = 1000, #Number of warmup samples to discard for each chain
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
out.

``` r

data("historical_posterior", package = "MBE")
head(historical_posterior)
#>    alphaCEonSur1Sur2 b1CEonSur1Sur2 b2CEonSur1Sur2 SigSqCEonSur1Sur2
#>                <num>          <num>          <num>             <num>
#> 1:      -0.060657298     -0.2923658    -0.02577700      0.0028743436
#> 2:      -0.050034948     -0.3666507    -0.02581211      0.0056530292
#> 3:      -0.021949313     -0.3229995    -0.02961769      0.0057219772
#> 4:      -0.003291189     -0.3301485    -0.02263008      0.0004533587
#> 5:      -0.045876147     -0.3155187    -0.02901523      0.0026134902
#> 6:      -0.023826164     -0.3730042    -0.02878529      0.0020488892
#>    alphaSur1onSur2  bSur1onSur2 SigSqSur1onSur2    muSur2 sigSqSur2
#>              <num>        <num>           <num>     <num>     <num>
#> 1:       0.4020287 -0.015995178       0.5994535 -1.708699  51.72702
#> 2:       0.3225021  0.002805734       0.5720484 -2.771109  47.25131
#> 3:       0.4069362  0.022673745       0.5007477 -2.547149  66.52862
#> 4:       0.3851224  0.022911445       0.5374545 -1.564569  53.40327
#> 5:       0.4598896  0.019820946       0.3365655 -1.126421  47.99716
#> 6:       0.4878410  0.029780083       0.3707338 -1.543360  54.41675
```

The column names of this posterior datasets are different from what’s
used in the manuscript. Here’s the mapping:

| Manuscript name         | Programming name    |
|:------------------------|:--------------------|
| $`\beta_0`$             | `alphaCEonSur1Sur2` |
| $`\beta_1`$             | `b1CEonSur1Sur2`    |
| $`\beta_2`$             | `b2CEonSur1Sur2`    |
| $`\lambda^2_\theta`$    | `SigSqCEonSur1Sur2` |
| $`\alpha_0`$            | `alphaSur1onSur2`   |
| $`\alpha_1`$            | `bSur1onSur2`       |
| $`\lambda^2_\gamma`$    | `SigSqSur1onSur2`   |
| $`\mu_{\gamma 2}`$      | `muSur2`            |
| $`\sigma^2_{\gamma 2}`$ | `sigSqSur2`         |

This posterior distribution is used as prior for estimating the
treatment effect on the clinical endpoint for the held-out trial. The
held-out trial is converted to the naming convention used by
[`MBE()`](https://hyejung0.github.io/MBE/reference/MBE.md):

``` r

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
posterior mixture. These quantities are calculated using weighted
sampling method shown in Appendix B of the manuscript. `post_psi0`
contains resampled posterior draws in the fixed order clinical endpoint,
surrogate 1, and surrogate 2. `weight_data$norm_w` is the set of
normalized weights, $`\{ \tilde{w}^{(b)} \}`$. `importance_diagnostics`
helps identify whether those weights are concentrated in only a small
number of historical draws.

## Historical-model fit assessment

Fitting and assessment are separate steps. Choose the historical-model
intercept and uncertainty prior, then fit the model to all but one trial
using
[`historical_model_fit_2surrogates()`](https://hyejung0.github.io/MBE/reference/historical_model_fit_2surrogates.md).
Repeat for each held-out trial. The two assessment functions,
[`loo_cv_rmse()`](https://hyejung0.github.io/MBE/reference/loo_cv_rmse.md)
and
[`tolerance_interval_coverage()`](https://hyejung0.github.io/MBE/reference/tolerance_interval_coverage.md),
only summarize inputs supplied by the user; they never fit models or
select priors.

This is internal leave-one-trial-out assessment: each held-out trial
serves as a pseudo-new trial. The package uses simulated trials, so its
results are not expected to match the manuscript’s CKD trial results.

### Prepare the fits and point estimates yourself

The following example shows the fitting loop explicitly. The settings
below are illustrative and can be changed independently of the
assessment functions. `niter` is the number of saved post-warmup draws
per chain. Run this expensive step only when new fits are needed; it is
not run when building the vignette.

``` r

historical_fields <- c(
  "CE_est", "CE_se", "Sur1_est", "Sur1_se", "Sur2_est", "Sur2_se",
  "Cor_CE_Sur1", "Cor_CE_Sur2", "Cor_Sur1_Sur2"
)
trials <- as.data.frame(trial_sim_dat)
historical_data <- trials[, historical_fields]

# Your choices for this analysis; the metric functions do not set these.
use_random_intercept <- TRUE
uncertainty_prior <- "inverse_gamma"
use_diffuse_surrogate_prior <- FALSE
center_clinical_intercept <- FALSE

paired_effects <- data.frame(
  observed = trials$CE_est,
  estimated = rep(NA_real_, nrow(trials)),
  row.names = as.character(trials$trial_id)
)
loo_draws <- vector("list", nrow(trials))
names(loo_draws) <- as.character(trials$trial_id)
loo_diagnostics <- vector("list", nrow(trials))

for (i in seq_len(nrow(trials))) {
  historical_fit <- historical_model_fit_2surrogates(
    data = historical_data[-i, , drop = FALSE],
    random_intercept = use_random_intercept,
    prior_for_uncertainty = uncertainty_prior,
    nchains = 4, ncores = 4, niter = 2000, nwarmup = 1000,
    seed = 2026 + i - 1,
    adapt_delta = 0.95, max_treedepth = 15
  )
  loo_draws[[i]] <- historical_fit$posterior_samples
  loo_diagnostics[[i]] <- list(
    sampler = historical_fit$diagnostics,
    convergence = historical_fit$fit$summary()
  )

  held_out <- list(
    ClnEst = trials$CE_est[i], ClnSE = trials$CE_se[i],
    Sur1Est = trials$Sur1_est[i], Sur1SE = trials$Sur1_se[i],
    Sur2Est = trials$Sur2_est[i], Sur2SE = trials$Sur2_se[i],
    R1Clin = trials$Cor_CE_Sur1[i],
    R2Clin = trials$Cor_CE_Sur2[i], R12 = trials$Cor_Sur1_Sur2[i]
  )
  set.seed(2026 + i - 1)
  updated <- MBE(
    mcmc_dat = loo_draws[[i]], sample_dat = held_out,
    diffuse = use_diffuse_surrogate_prior,
    intercept0 = center_clinical_intercept
  )
  paired_effects$estimated[i] <- updated$post_mean[1, 1]
}

# Inspect diagnostics and resolve fitting issues before interpreting the metric.
loo_diagnostics[[1]]
loo_cv_rmse(paired_effects)
```

If you already have fits, reuse their posterior draws in the
[`MBE()`](https://hyejung0.github.io/MBE/reference/MBE.md) calls; there
is no need to repeat the historical fitting step. Saved records from the
earlier development workflow can still be read with
`readRDS(path)$posterior_draws`.

For example, if the original 66 saved records are on your computer, load
them directly instead of running the fitting loop:

``` r

trials <- as.data.frame(trial_sim_dat)
loo_draws <- lapply(as.character(trials$trial_id), function(id) {
  record <- readRDS(file.path(
    "data-raw/loo-random-inverse-gamma", paste0("heldout_trial_", id, ".rds")
  ))
  stopifnot(identical(as.character(record$held_out_trial), id))
  record$posterior_draws
})
names(loo_draws) <- as.character(trials$trial_id)
use_diffuse_surrogate_prior <- FALSE
center_clinical_intercept <- FALSE
```

### Calculate RMSE from a two-column table

Supply one row per trial and two columns named `observed` and
`estimated`. A data frame or numeric matrix is accepted, with columns in
either order. Match trial IDs before assembling the pairs. The estimate
is chosen by the user; the example uses the posterior mean from
[`MBE()`](https://hyejung0.github.io/MBE/reference/MBE.md), an estimate
of the true effect rather than a known true effect.

``` r

data("loo_assessment_by_trial", package = "MBE")
paired_effects <- data.frame(
  observed = loo_assessment_by_trial$observed_clinical,
  estimated = loo_assessment_by_trial$posterior_mean_clinical
)
loo_cv_rmse(paired_effects)
#> [1] 0.2385382
loo_cv_rmse(as.matrix(paired_effects))
#> [1] 0.2385382
```

The calculation is `sqrt(mean((observed - estimated)^2))`. With a full
MBE update, the held-out clinical estimate contributes to its own
posterior mean. This measures agreement after updating, rather than
prediction from surrogates alone. Calling the metric LOO-CV assumes the
corresponding trial was excluded from historical fitting; the metric
function cannot verify that step.

### Generate the clinical predictive distribution and check capture

[`clinical_predictive_distribution()`](https://hyejung0.github.io/MBE/reference/clinical_predictive_distribution.md)
takes historical posterior draws and the held-out trial’s surrogate
estimates, standard errors, and sampling correlations. It generates
draws for the observed clinical estimate conditional on the two
surrogate estimates. The clinical estimate itself is ignored even if
`ClnEst` is present in the input list. This example reuses
`historical_posterior` and `new_trial` from above; trial 1 was excluded
from that historical fit:

``` r

predictive_draws <- clinical_predictive_distribution(
  mcmc_dat = historical_posterior, sample_dat = new_trial,
  diffuse = FALSE, intercept0 = FALSE,
  n_draws = 4000, seed = 2026
)
stats::quantile(predictive_draws, c(0.025, 0.975))
#>       2.5%      97.5% 
#> -0.3745658  0.1562384
```

[`tolerance_interval_coverage()`](https://hyejung0.github.io/MBE/reference/tolerance_interval_coverage.md)
calls this prediction function internally. Supply the desired lower and
upper quantile probabilities: `c(0.025, 0.975)` for a 95% interval, or
one row per interval in a two-column matrix. Asymmetric intervals are
also supported.

``` r

interval_probs <- rbind(
  c(0.025, 0.975),  # 95%
  c(0.05, 0.95),    # 90%
  c(0.075, 0.925)   # 85%
)
tolerance_interval_coverage(
  observed = new_trial$ClnEst,
  mcmc_dat = historical_posterior, sample_dat = new_trial,
  probs = interval_probs,
  diffuse = FALSE, intercept0 = FALSE,
  n_draws = 4000, seed = 2026
)
#>   lower_prob upper_prob level   observed      lower      upper covered
#> 1      0.025      0.975    95 -0.1337531 -0.3745658 0.15623836    TRUE
#> 2      0.050      0.950    90 -0.1337531 -0.3326389 0.10673561    TRUE
#> 3      0.075      0.925    85 -0.1337531 -0.3078896 0.07480282    TRUE
```

The result has one row per probability pair: `lower_prob`, `upper_prob`,
`level` (as a percentage), `observed`, `lower`, `upper`, and logical
`covered`. Bounds use `stats::quantile(type = 7)`, and both endpoints
count as covered. The function generates predictive draws once per
trial, then uses the same draws for all requested intervals. A fixed
seed makes simulation reproducible.

For multiple trials, reuse the `loo_draws` list from the fitting
example, match it to the trial IDs, and combine single-trial results:

``` r

stopifnot(
  !is.null(names(loo_draws)), !anyDuplicated(names(loo_draws)),
  !anyDuplicated(trials$trial_id),
  setequal(as.character(trials$trial_id), names(loo_draws))
)
loo_draws <- loo_draws[as.character(trials$trial_id)]
capture <- do.call(rbind, lapply(seq_len(nrow(trials)), function(i) {
  held_out <- list(
    ClnSE = trials$CE_se[i],
    Sur1Est = trials$Sur1_est[i], Sur1SE = trials$Sur1_se[i],
    Sur2Est = trials$Sur2_est[i], Sur2SE = trials$Sur2_se[i],
    R1Clin = trials$Cor_CE_Sur1[i],
    R2Clin = trials$Cor_CE_Sur2[i], R12 = trials$Cor_Sur1_Sur2[i]
  )
  result <- tolerance_interval_coverage(
    observed = trials$CE_est[i],
    mcmc_dat = loo_draws[[i]], sample_dat = held_out,
    probs = interval_probs, n_draws = 8000, seed = 2026 + i - 1,
    diffuse = use_diffuse_surrogate_prior,
    intercept0 = center_clinical_intercept
  )
  data.frame(trial_id = trials$trial_id[i], result)
}))
stats::aggregate(covered ~ lower_prob + upper_prob, data = capture, FUN = mean)
```

Prediction conditions the joint observed-effect distribution on the two
surrogate estimates. It incorporates residual clinical heterogeneity,
clinical sampling variance, surrogate uncertainty, and the supplied
sampling correlations. Historical draws are weighted by the likelihood
of the observed surrogates. This is different from passing latent-effect
draws from `MBE()$post_psi0`, which already condition on the observed
clinical estimate.

### Archived example results

The completed 66-trial assessment remains available as two compact
datasets. Their provenance is the random-intercept, inverse-gamma
specification; neither new metric requires that specification. The
archived coverage used an empirical percentile rule and an unweighted
surrogate-update simulation. The new function uses joint conditioning,
surrogate-based weights, and inclusive quantile bounds. Those
differences can change coverage. Archived results are preserved, not
recalculated or relabeled as results of the new prediction function.

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
