# Calculate LOO-CV RMSE and Tolerance-Interval Coverage

Calculates the historical-model assessment measures described in the MBE
manuscript from posterior draws obtained after holding out each trial.
LOO-CV RMSE uses the posterior clinical mean from a full-carryover MBE
update. Tolerance coverage predicts the held-out observed clinical
estimate using its two observed surrogate estimates without conditioning
on the observed clinical estimate itself.

## Usage

``` r
loo_cv_model_assessment(loo_posteriors, data, trial_id = NULL, seed = NULL)
```

## Arguments

- loo_posteriors:

  A list with one posterior-draw data frame or matrix per held-out
  trial, or a character vector of `.rds` paths created by
  [`fit_loo_historical_models()`](https://hyejung0.github.io/MBE/reference/fit_loo_historical_models.md).
  The output directory itself may also be supplied, in which case its
  saved trial files are discovered and matched using their metadata. A
  single data frame or matrix is accepted when `data` contains one
  trial. Named inputs are matched to trial IDs; unnamed inputs must
  follow the row order of `data`.

- data:

  A data frame or matrix containing the held-out trial observations. It
  must include the nine historical-model fields; additional columns are
  allowed.

- trial_id:

  Optional name of the unique trial-ID column. When `NULL`, `trial_id`
  is used if present; otherwise row numbers are used.

- seed:

  Optional positive integer used for reproducible posterior predictive
  simulation. The caller's random-number state is preserved.

## Value

A list with two data frames:

- per_trial:

  Held-out observations, posterior clinical means, squared errors,
  predictive percentiles and intervals, and coverage indicators.

- summary:

  Number of trials, LOO-CV RMSE, and aggregate 90% and 95%
  tolerance-interval coverage.

## Details

For each MCMC draw, the two true surrogate effects are updated using the
two held-out surrogate estimates. A predictive draw for the observed
clinical estimate is then generated with variance
`SigSqCEonSur1Sur2 + ClnSE^2`. This corrects the original analysis
script, which inadvertently omitted `SigSqCEonSur1Sur2` from the
predictive variance.

The 90% and 95% coverage indicators follow the manuscript calculation:
the empirical percentile of the observed clinical estimate must fall
within 0.05–0.95 or 0.025–0.975, respectively.

## Examples

``` r
data("historical_posterior")
data("trial_sim_dat")

# The packaged historical posterior was fit with trial 1 held out, so it can
# demonstrate the calculation for one trial. Aggregate coverage requires one
# separately fitted posterior for every held-out trial.
one_trial <- loo_cv_model_assessment(
  loo_posteriors = list(`1` = historical_posterior[1:100, ]),
  data = trial_sim_dat[1, ],
  seed = 1
)
one_trial$per_trial
#>   trial_id observed_clinical posterior_mean_clinical squared_error
#> 1        1        -0.1337531              -0.1223965   0.000128971
#>   predictive_percentile predictive_lower_95 predictive_lower_90
#> 1                  0.47          -0.4625712          -0.4314274
#>   predictive_upper_90 predictive_upper_95 covered_95 covered_90
#> 1           0.1749716            0.213356       TRUE       TRUE
one_trial$summary
#>   n_trials loo_cv_rmse coverage_95 coverage_90
#> 1        1  0.01135654           1           1
```
