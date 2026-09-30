# Fit Random-Intercept Inverse-Gamma Models for Leave-One-Out Assessment

Fits the package's initial manuscript specification once for each
held-out trial: a random clinical intercept with inverse-gamma priors
for variance parameters. Each fit uses all remaining trials. Only the
nine posterior columns needed for MBE assessment are saved, rather than
complete CmdStan fit objects or chain CSV files.

## Usage

``` r
fit_loo_historical_models(
  data,
  output_dir,
  trial_id = NULL,
  seed = 1L,
  overwrite = FALSE,
  nchains = 4L,
  ncores = 1L,
  niter = 2000L,
  nwarmup = 1000L,
  show_messages = TRUE,
  ...
)
```

## Arguments

- data:

  A data frame or matrix containing the nine historical-model fields
  documented in
  [`historical_model_fit_2surrogates()`](https://hyejung0.github.io/MBE/reference/historical_model_fit_2surrogates.md).
  Additional columns, such as `trial_id` and simulation truths, are
  allowed.

- output_dir:

  Directory in which to save one compact `.rds` file per held-out trial.
  A directory under `data-raw/` is recommended; the full files should
  not be included in the installed package.

- trial_id:

  Optional name of a column containing unique trial IDs. When `NULL`,
  `trial_id` is used if present; otherwise row numbers are used.

- seed:

  Positive integer base seed. Trial `i` uses `seed + i - 1`.

- overwrite:

  Logical; overwrite existing trial files when `TRUE`.

- nchains:

  Number of MCMC chains for each historical fit.

- ncores:

  Number of chains to run in parallel within each fit.

- niter:

  Number of post-warmup draws per chain.

- nwarmup:

  Number of warmup iterations per chain.

- show_messages:

  Logical; show CmdStan sampling messages.

- ...:

  Additional arguments forwarded to
  [`historical_model_fit_2surrogates()`](https://hyejung0.github.io/MBE/reference/historical_model_fit_2surrogates.md)
  and then to `CmdStanModel$sample()`, such as `adapt_delta`,
  `max_treedepth`, and `refresh`.

## Value

Invisibly, a named character vector of saved `.rds` paths. Names are the
held-out trial identifiers. Existing files are returned without
refitting when `overwrite = FALSE`, after their trial identifier, seed,
model specification, and required posterior columns are validated.

## Examples

``` r
if (FALSE) { # \dontrun{
paths <- fit_loo_historical_models(
  trial_sim_dat,
  output_dir = "data-raw/loo-random-inverse-gamma",
  seed = 2026,
  nchains = 4,
  ncores = 4,
  niter = 2000,
  nwarmup = 1000
)
} # }
```
