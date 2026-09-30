# Estimate an MBE Posterior with Two Surrogate Endpoints

Uses new-trial estimates and historical posterior draws to estimate the
posterior distribution of an MBE by importance sampling. The model is
fixed to one clinical endpoint and exactly two surrogate endpoints.

## Usage

``` r
MBE(mcmc_dat, sample_dat, diffuse = TRUE, diffuse_se = 100, intercept0 = TRUE)
```

## Arguments

- mcmc_dat:

  A data frame or matrix containing posterior draws from the
  two-surrogate historical model.

- sample_dat:

  A named list containing exactly nine scalar values: `ClnEst`, `ClnSE`,
  `Sur1Est`, `Sur1SE`, `Sur2Est`, `Sur2SE`, `R1Clin`, `R2Clin`, and
  `R12`.

- diffuse:

  A logical value indicating whether to use a diffuse prior. Defaults to
  TRUE.

- diffuse_se:

  A positive numeric value specifying the standard deviation of the
  diffuse prior for the two surrogate effects.

- intercept0:

  A logical value indicating whether to center the intercept term of the
  meta-regression to zero. Defaults to TRUE.

## Value

A list containing posterior means, posterior covariance, weighted
quantiles for the clinical and two surrogate effects, posterior draws,
the corresponding importance weights, and importance-sampling
diagnostics. The diagnostic vector reports the effective sample size,
effective sample size relative to the number of historical draws, and
maximum normalized weight. Endpoint components are always ordered as
clinical endpoint, surrogate 1, and surrogate 2.

## Examples

``` r
# Use the example data provided in the package.
# Load the historical posterior samples where the model was fit
# with random intercept and inverse gamma prior for variance parameters
# on the 65/66 simulated trials in `trial_sim_dat` dataset.
# It's the first trial that was left out as if it was a new trial.
# We demonstrate estimating the posterior distribution of MBE on this first
# row that was left out.
data("historical_posterior")
data("trial_sim_dat")
one_sim_dat<-list(

# treatment effect on CE
ClnEst=trial_sim_dat$CE_est[1],
ClnSE=trial_sim_dat$CE_se[1],

# treatment effect on chronic slope
Sur1Est=trial_sim_dat$Sur1_est[1],
Sur1SE=trial_sim_dat$Sur1_se[1],

# treatment effect on acute slope
Sur2Est=trial_sim_dat$Sur2_est[1],
Sur2SE=trial_sim_dat$Sur2_se[1],

# correlation between CE and chronic slope
R1Clin=trial_sim_dat$Cor_CE_Sur1[1],

# correlation between CE and acute slope
R2Clin=trial_sim_dat$Cor_CE_Sur2[1],

# correlation between chronic slope and acute slope
R12=trial_sim_dat$Cor_Sur1_Sur2[1]
)

MBE_distribution<-MBE(
mcmc_dat = historical_posterior,
sample_dat=one_sim_dat,
diffuse_se = 100,
diffuse = TRUE,
intercept0 = TRUE
)

head(MBE_distribution$post_psi0)
#>         clinical surrogate1 surrogate2
#> [1,] -0.01742418  0.4358344  -7.647448
#> [2,] -0.05070114  0.3625249  -5.290546
#> [3,] -0.21549799  0.9284195  -4.653444
#> [4,] -0.22270912  1.0586468  -4.846612
#> [5,] -0.02885444  0.6432777  -8.217527
#> [6,] -0.22661395  0.9851503  -4.977536
```
