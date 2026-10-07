# Trial-Level Leave-One-Out Historical-Model Assessment

Trial-level results underlying `loo_assessment_summary`. The assessment
uses the random-intercept historical model with inverse-gamma variance
priors and seed 2026.

## Usage

``` r
loo_assessment_by_trial
```

## Format

A data frame with 66 rows and 11 variables:

- trial_id:

  Character. Identifier of the held-out trial.

- observed_clinical:

  Numeric. Observed clinical treatment-effect estimate.

- posterior_mean_clinical:

  Numeric. Posterior mean for the held-out clinical treatment effect.

- squared_error:

  Numeric. Squared error of the posterior mean.

- predictive_percentile:

  Numeric. Empirical predictive percentile of the observed clinical
  estimate.

- predictive_lower_95:

  Numeric. Lower endpoint of the 95 percent predictive interval.

- predictive_lower_90:

  Numeric. Lower endpoint of the 90 percent predictive interval.

- predictive_upper_90:

  Numeric. Upper endpoint of the 90 percent predictive interval.

- predictive_upper_95:

  Numeric. Upper endpoint of the 95 percent predictive interval.

- covered_95:

  Logical. Whether the 95 percent interval covered the observed clinical
  estimate.

- covered_90:

  Logical. Whether the 90 percent interval covered the observed clinical
  estimate.

## Source

Archived 66-fit analysis of `trial_sim_dat` with seed 2026, saved in
`data-raw/loo-random-inverse-gamma/loo_assessment.rds`. See the
preparation script `data-raw/loo_assessment.R` in the source repository.

## Details

These archived results use the earlier empirical-percentile coverage
calculation; see
[loo_assessment_summary](https://hyejung0.github.io/MBE/reference/loo_assessment_summary.md)
for provenance and its differences from the current prediction and
coverage functions. The paired effects can still be supplied to
[`loo_cv_rmse()`](https://hyejung0.github.io/MBE/reference/loo_cv_rmse.md)
without refitting any model.
