# Changelog

## MBE 0.0.0.9000

- Separate historical fitting from assessment. Use
  [`historical_model_fit_2surrogates()`](https://hyejung0.github.io/MBE/reference/historical_model_fit_2surrogates.md)
  directly on each training set, then
  [`loo_cv_rmse()`](https://hyejung0.github.io/MBE/reference/loo_cv_rmse.md)
  on an observed/estimated table and
  [`tolerance_interval_coverage()`](https://hyejung0.github.io/MBE/reference/tolerance_interval_coverage.md)
  on existing posterior draws and held-out sampling inputs.
- Remove the development functions `fit_loo_historical_models()` and
  `loo_cv_model_assessment()`. Assessment no longer chooses an intercept
  or variance prior, fits historical models, or updates MBE posteriors.
- Add
  [`clinical_predictive_distribution()`](https://hyejung0.github.io/MBE/reference/clinical_predictive_distribution.md),
  called internally by coverage. It conditions the joint observed-effect
  model on surrogate estimates and weights historical draws by surrogate
  likelihood, accounting for clinical- surrogate sampling correlations
  without using the observed clinical estimate.
- Coverage accepts a pair of quantile probabilities (`c(0.025, 0.975)`)
  or a two-column matrix of pairs. It tests inclusive quantile bounds.
  Previous assessment used empirical percentiles and an unweighted
  surrogate-update simulation; new coverage need not match archived
  results.
- Preserve the completed 66-trial assessment tables and locally saved
  fits as archived examples; no refitting is required to use the new
  RMSE function.
