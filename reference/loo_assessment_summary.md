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

Archived 66-fit analysis of `trial_sim_dat` with seed 2026, saved in
`data-raw/loo-random-inverse-gamma/loo_assessment.rds`. See the
preparation script `data-raw/loo_assessment.R` in the source repository.

## Details

Archived results from the earlier development assessment. RMSE uses
full-carryover MBE posterior means, incorporating the held-out clinical
estimate. Coverage used empirical-percentile indicators and an
unweighted surrogate-update simulation. The current
[`clinical_predictive_distribution()`](https://hyejung0.github.io/MBE/reference/clinical_predictive_distribution.md)
uses joint conditioning and surrogate-likelihood weights, and
[`tolerance_interval_coverage()`](https://hyejung0.github.io/MBE/reference/tolerance_interval_coverage.md)
checks inclusive quantile bounds. This table is retained unchanged and
is not a result of the current coverage function. These are internal
cross-validation results on simulated data, not external validation.
Model specification here describes only this archived dataset; the
current assessment functions do not prescribe an intercept or prior.
