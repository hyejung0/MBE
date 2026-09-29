# Simulated Trial-Level Data for Two Surrogate Endpoints

Simulated summary data for 66 randomized trials. Each trial has one
clinical endpoint and exactly two surrogate endpoints (chronic and acute
eGFR slopes), along with their standard errors, correlations, and
simulation truths.

## Usage

``` r
trial_sim_dat
```

## Format

A data frame with 66 rows (one row per trial) and 13 variables:

- trial_id:

  Integer. Index for trial.

- CE_est:

  Numeric. Observed estimated treatment effect on the clinical endpoint.

- Sur1_est:

  Numeric. Observed estimated treatment effect on the first surrogate
  endpoint.

- Sur2_est:

  Numeric. Observed estimated treatment effect on the second surrogate
  endpoint.

- CE_se:

  Numeric. Standard error of the observed estimated treatment effect on
  the clinical endpoint.

- Sur1_se:

  Numeric. Standard error of the observed estimated treatment effect on
  the first surrogate endpoint.

- Sur2_se:

  Numeric. Standard error of the observed estimated treatment effect on
  the second surrogate endpoint.

- Cor_CE_Sur1:

  Numeric. Estimated correlation between clinical endpoint and 1st
  surrogate endpoint.

- Cor_CE_Sur2:

  Numeric. Estimated correlation between clinical endpoint and 2nd
  surrogate endpoint.

- Cor_Sur1_Sur2:

  Numeric. Estimated correlation between 1st surrogate endpoint and 2nd
  surrogate endpoint.

- theta_CE:

  Numeric. True treatment effect on the clinical endpoint.

- theta_Sur1:

  Numeric. True treatment effect on the 1st surrogate endpoint.

- theta_Sur2:

  Numeric. True treatment effect on the 2nd surrogate endpoint.

## Source

Simulated by Yizhen Xu

## Author

Yizhen Xu <yizhen.xu@utah.edu>
