# Aggregate Leave-One-Out Historical-Model Assessment

A one-row summary of leave-one-out cross-validation for the
random-intercept historical model with inverse-gamma variance priors.
Each of the 66 simulated trials was held out once. The
tolerance-interval calculation includes both residual clinical
heterogeneity and the held-out trial's clinical sampling variance.

## Usage

``` r
loo_assessment_summary
```

## Format

A data frame with 1 row and 4 variables:

- n_trials:

  Integer. Number of held-out trials.

- loo_cv_rmse:

  Numeric. Root mean squared error between held-out clinical estimates
  and their posterior means.

- coverage_95:

  Numeric. Proportion of held-out clinical estimates covered by the 95
  percent predictive interval.

- coverage_90:

  Numeric. Proportion of held-out clinical estimates covered by the 90
  percent predictive interval.

## Source

Derived from `trial_sim_dat` using
[`fit_loo_historical_models()`](https://hyejung0.github.io/MBE/reference/fit_loo_historical_models.md)
and
[`loo_cv_model_assessment()`](https://hyejung0.github.io/MBE/reference/loo_cv_model_assessment.md)
with seed 2026.
