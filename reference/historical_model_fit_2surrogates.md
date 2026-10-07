# Fit the Historical Model for Exactly Two Surrogate Endpoints

Fits a trial-level meta-analysis using a two-stage random-effects
Bayesian model for exactly two surrogate endpoints and one definitive
clinical endpoint, as described in Lee et al. (2026).

## Usage

``` r
historical_model_fit_2surrogates(
  data,
  random_intercept = TRUE,
  prior_for_uncertainty = "half_normal",
  output_dir = tempdir(),
  show_messages = TRUE,
  nchains = 4,
  ncores = 1,
  niter = 2000,
  nwarmup = 1000,
  ...
)
```

## Arguments

- data:

  A data frame or matrix with one row per trial and exactly these nine
  numeric columns (order does not matter): `CE_est`, `CE_se`,
  `Sur1_est`, `Sur1_se`, `Sur2_est`, `Sur2_se`, `Cor_CE_Sur1`,
  `Cor_CE_Sur2`, and `Cor_Sur1_Sur2`. Data for a third surrogate
  endpoint are not accepted.

- random_intercept:

  logical indicating whether to include intercept in regression
  modeling.

- prior_for_uncertainty:

  character indicating the prior distribution for the uncertainty
  parameters. Options are "half_normal" for half-normal prior on the
  standard deviation scale and "inverse_gamma" for inverse-gamma prior
  on the variance scale. Default is "half_normal".

- output_dir:

  Character string. The folder path where CmdStan should save the raw
  MCMC chain CSV files. Defaults to
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html), which saves files
  to a temporary session folder that is automatically deleted when R
  closes. To keep the CSV files permanently, provide a local directory
  path (e.g., `"./"` for the current working directory, or
  `"./mcmc_output"`).

- show_messages:

  Logical; whether to print progress to console. Default is TRUE.

- nchains:

  number of chains for MCMC sampling. Default to 4.

- ncores:

  number of chains to run in parallel. Default to 1 (sequential
  processing). If set to a value greater than 1, the function will use
  parallel processing to run multiple chains simultaneously.

- niter:

  number of post-warmup samples to save for each chain. Default is 2000.

- nwarmup:

  number of warmup iterations for MCMC sampling.

- ...:

  Additional arguments forwarded to `cmdstanr`'s `$sample()` method
  (e.g., `adapt_delta`, `max_treedepth`, `seed`, `refresh`). The
  `cmdstanr`'s default values are used for `adapt_delta = 0.9` and
  `max_treedepth = 12`. In addition, `chains`, `parallel_chains`,
  `iter_sampling`, and `iter_warmup` are set to the values of `nchains`,
  `ncores`, `niter`, and `nwarmup`, respectively, unless overridden in
  `...`.

## Value

A list containing the fitted model object and other relevant
information:

- diagnostics:

  Diagnostic summary table from `fit$diagnostic_summary()`.

- posterior_samples:

  Sampled posterior distribution draws as a data frame.

- summary:

  A summary of the posterior distribution of parameters of interest.

- loo:

  Efficient approximate leave-one-out (LOO) cross-validation from
  [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html).

- waic:

  Widely applicable information criterion (WAIC) calculated with
  [`loo::waic()`](https://mc-stan.org/loo/reference/waic.html).

- fit:

  The fitted `CmdStanMCMC` object from `cmdstanr`.

## Details

This function fits a trial-level meta-analysis model using a 2-stage
random effects Bayesian model under chronic kidney disease (CKD) context
first introduced by Lee et al. (2026). For a true treatment effect on
clinical endpoint \\\theta_i\\, and true treatment effect on two
surrogate endpoints \\\gamma\_{i, 1}\\ and \\\gamma\_{i, 2}\\ for the
\\i\\-th trial, the model is specified as follows: \$\$\hat{\psi}\_i
\mid \psi_i \sim N_3(\psi_i, \Sigma\_{y, i})\$\$ where \\\psi_i =
(\theta_i, \gamma\_{i, 1}, \gamma\_{i, 2})^T\\ is the true treatment
effect vector for the \\i\\-th trial, \\\hat{\psi}\_i\\ is the observed
(estimated) treatment effect vector for the \\i\\-th trial, and
\\\Sigma\_{y, i}\\ is the covariance matrix of the corresponding
\\\hat{\psi}\_i\\.

The density of \\\psi\\ can be expressed as a sequence of conditional
densities: \$\$\gamma\_{i, 2} \sim N(\mu_2, \sigma_2^2)\$\$
\$\$\gamma\_{i, 1} \mid \gamma\_{i, 2} \sim N(\alpha\_{\gamma_1} +
\beta\_{\gamma_1} \cdot \gamma\_{i, 2}, \lambda^2\_{\gamma_1})\$\$
\$\$\theta_i \mid \gamma\_{i, 1}, \gamma\_{i, 2} \sim N(\alpha\_\theta +
\beta\_{\gamma_1} \cdot \gamma\_{i, 1} + \beta\_{\gamma_2} \cdot
\gamma\_{i, 2}, \lambda^2\_\theta)\$\$

This function fits the model using MCMC sampling on historical RCTs and
returns the posterior distribution of the model parameters.

The `loo` and `waic` components are conditional model-fit measures based
on the manuscript's Stan log-likelihood. They are distinct from the
explicit internal leave-one-trial-out assessment. For that assessment,
call this fitting function separately for each set of training trials,
omitting the held-out trial. Then use
[`loo_cv_rmse()`](https://hyejung0.github.io/MBE/reference/loo_cv_rmse.md)
on paired effects and
[`tolerance_interval_coverage()`](https://hyejung0.github.io/MBE/reference/tolerance_interval_coverage.md)
with the already fitted posterior draws and held-out sampling inputs.
The latter calls
[`clinical_predictive_distribution()`](https://hyejung0.github.io/MBE/reference/clinical_predictive_distribution.md).
The assessment functions do not fit models or choose priors.

## Examples

``` r
if (FALSE) { # \dontrun{
data("trial_sim_dat", package = "MBE")

historical_data <- trial_sim_dat[, c(
  "CE_est", "CE_se", "Sur1_est", "Sur1_se", "Sur2_est", "Sur2_se",
  "Cor_CE_Sur1", "Cor_CE_Sur2", "Cor_Sur1_Sur2"
)]

fit <- historical_model_fit_2surrogates(
  data = historical_data,
  random_intercept = TRUE,
  prior_for_uncertainty = "half_normal",
  nchains = 4,
  ncores = 4,
  niter = 2000,
  nwarmup = 1000,
  adapt_delta = 0.95,  # Target acceptance rate (resolves divergences)
  max_treedepth = 15,  # Tree depth limit (resolves treedepth saturation)
  refresh = 250        # Print updates every 250 iterations
)
} # }
```
