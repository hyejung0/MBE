# Interim Analysis Example

`interim_sim_dat` contains simulated interim analyses for eight
scenarios. Each row has one clinical endpoint and exactly two surrogate
endpoints: chronic (`Sur1`) and acute (`Sur2`) eGFR slopes. Some early
analyses have missing clinical estimates because too few clinical events
had occurred; those rows are not yet eligible for an MBE update.

``` r

new_trial_fields <- c(
  "ClnEst", "ClnSE",
  "Sur1Est", "Sur1SE",
  "Sur2Est", "Sur2SE",
  "R1Clin", "R2Clin", "R12"
)

interim_data <- as.data.frame(interim_sim_dat)
eligible <- stats::complete.cases(interim_data[, new_trial_fields])
table(eligible)
#> eligible
#> FALSE  TRUE 
#>   123   405
```

The first complete interim analysis can be updated using the historical
posterior supplied with the package:

``` r

first_complete <- which(eligible)[1]
new_trial <- as.list(interim_data[first_complete, new_trial_fields])

set.seed(1)
interim_mbe <- MBE(
  mcmc_dat = historical_posterior,
  sample_dat = new_trial,
  diffuse = TRUE,
  diffuse_se = 100,
  intercept0 = TRUE
)

interim_data[first_complete, c("case", "analysis.month")]
#>    case analysis.month
#> 19    1             25
interim_mbe$post_mean
#>                   [,1]
#> clinical    0.00675136
#> surrogate1 -0.00804139
#> surrogate2 -0.14778649
interim_mbe$post_quantiles_clinical
#> quantile_0.025  quantile_0.05   quantile_0.1  quantile_0.25   quantile_0.5 
#>   -0.304055269   -0.252149305   -0.195916295   -0.100410823    0.007711854 
#>  quantile_0.75   quantile_0.9  quantile_0.95 quantile_0.975 
#>    0.113076436    0.210184273    0.270522806    0.329201746
```

For a sequence of interim analyses, repeat the construction of
`new_trial` for each complete row. The historical posterior can be
reused, but the treatment effect estimates, standard errors, and
correlations must all correspond to the same interim time point.
