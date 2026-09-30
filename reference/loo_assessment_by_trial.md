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

Derived from `trial_sim_dat` using
[`fit_loo_historical_models()`](https://hyejung0.github.io/MBE/reference/fit_loo_historical_models.md)
and
[`loo_cv_model_assessment()`](https://hyejung0.github.io/MBE/reference/loo_cv_model_assessment.md)
with seed 2026.
