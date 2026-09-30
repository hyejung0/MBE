# Simulation Example

## Simulated trial data

`trial_sim_dat` contains 66 simulated trials with one clinical endpoint
and exactly two surrogate endpoints. Each row contains treatment-effect
estimates, standard errors, sampling correlations, and the corresponding
simulation truths.

``` r

head(trial_sim_dat)
#>    trial_id     CE_est   Sur1_est  Sur2_est     CE_se    Sur1_se   Sur2_se
#>       <int>      <num>      <num>     <num>     <num>      <num>     <num>
#> 1:        1 -0.1337531  0.6746282 -5.616994 0.1345275 0.24675735 2.0208922
#> 2:        2  0.3126354  0.9536621 -7.708984 0.1812155 0.37357333 2.9033496
#> 3:        3 -0.7966332  1.6220178  3.565406 0.1539254 0.26112209 2.0994112
#> 4:        4 -0.3433848  1.3337097  3.450051 0.4338270 0.43152790 4.8085586
#> 5:        5 -0.1624861 -0.4845420  3.818558 0.4265650 0.43311521 4.8014942
#> 6:        6 -0.5146880  0.5651871  5.729746 0.1156836 0.08444745 0.8259361
#>    Cor_CE_Sur1 Cor_CE_Sur2 Cor_Sur1_Sur2    theta_CE theta_Sur1 theta_Sur2
#>          <num>       <num>         <num>       <num>      <num>      <num>
#> 1:  -0.5432718  -0.2525332       -0.0045 -0.18937028  0.7879757  -6.431669
#> 2:  -0.5444478  -0.2921709        0.0206  0.13645441  0.4905083  -6.672995
#> 3:  -0.5455985  -0.2470125        0.0052 -0.31625150  0.8574721   2.082759
#> 4:  -0.5474905  -0.0928973       -0.2283 -0.54795430  1.4566659   1.017808
#> 5:  -0.5228790  -0.1340568       -0.2185 -0.06528294 -0.7109980   5.768145
#> 6:  -0.1667742  -0.1519148       -0.1580 -0.32536177  0.4891863   4.465707
```

The column definitions are available with
[`?trial_sim_dat`](https://hyejung0.github.io/MBE/reference/trial_sim_dat.md).

## Fit the historical model

For this example, trial 1 is treated as a new trial and the historical
model is fit to the remaining 65 trials. The call below uses the
manuscript’s initial random-intercept, inverse-gamma specification. It
is not run while building the vignette because it requires CmdStan and
several minutes of computation.

``` r

historical_fit <- historical_model_fit_2surrogates(
  data = trial_sim_dat[-1, ],
  random_intercept = TRUE,
  prior_for_uncertainty = "inverse_gamma",
  nchains = 4,
  ncores = 4,
  niter = 2000,
  nwarmup = 1000,
  seed = 2026,
  adapt_delta = 0.95,
  max_treedepth = 15
)

historical_fit$diagnostics
historical_fit$summary

historical_posterior <- as.data.frame(historical_fit$posterior_samples)[, c(
  "alphaCEonSur1Sur2",
  "b1CEonSur1Sur2",
  "b2CEonSur1Sur2",
  "SigSqCEonSur1Sur2",
  "alphaSur1onSur2",
  "bSur1onSur2",
  "SigSqSur1onSur2",
  "muSur2",
  "sigSqSur2"
)]
```

The package’s `historical_posterior` dataset contains compatible
precomputed draws from this specification with trial 1 held out.

## Update the MBE

Construct the new-trial input from the held-out row:

``` r

new_trial <- list(
  ClnEst = trial_sim_dat$CE_est[1],
  ClnSE = trial_sim_dat$CE_se[1],
  Sur1Est = trial_sim_dat$Sur1_est[1],
  Sur1SE = trial_sim_dat$Sur1_se[1],
  Sur2Est = trial_sim_dat$Sur2_est[1],
  Sur2SE = trial_sim_dat$Sur2_se[1],
  R1Clin = trial_sim_dat$Cor_CE_Sur1[1],
  R2Clin = trial_sim_dat$Cor_CE_Sur2[1],
  R12 = trial_sim_dat$Cor_Sur1_Sur2[1]
)
```

Then update the historical distribution with the new-trial estimates:

``` r

set.seed(1)
mbe_fit <- MBE(
  mcmc_dat = historical_posterior,
  sample_dat = new_trial,
  diffuse = TRUE,
  diffuse_se = 100,
  intercept0 = TRUE
)

mbe_fit$post_mean
#>                  [,1]
#> clinical   -0.1013294
#> surrogate1  0.6837417
#> surrogate2 -5.4760921
mbe_fit$post_var
#>               clinical   surrogate1   surrogate2
#> clinical    0.01196417 -0.019592518 -0.091338433
#> surrogate1 -0.01959252  0.060139164 -0.009146439
#> surrogate2 -0.09133843 -0.009146439  3.981621682
mbe_fit$post_quantiles_clinical
#> quantile_0.025  quantile_0.05   quantile_0.1  quantile_0.25   quantile_0.5 
#>    -0.31485836    -0.28123031    -0.24074706    -0.17402273    -0.09799387 
#>  quantile_0.75   quantile_0.9  quantile_0.95 quantile_0.975 
#>    -0.02681263     0.03707755     0.07905239     0.12341721
mbe_fit$importance_diagnostics
#>          effective_sample_size relative_effective_sample_size 
#>                   3.943574e+03                   9.858936e-01 
#>      maximum_normalized_weight 
#>                   3.262557e-04
```

The rows of `post_mean` and `post_var` are ordered as the clinical
endpoint, surrogate 1, and surrogate 2. The posterior is a mixture over
historical-model draws; `post_mean` and `post_var` summarize its first
two moments rather than asserting that the mixture itself is a single
multivariate normal distribution.

Posterior samples and normalized importance weights are also available:

``` r

head(mbe_fit$post_psi0)
#>         clinical surrogate1 surrogate2
#> [1,]  0.08370764  0.3990721  -7.852443
#> [2,] -0.11376674  0.8685964  -5.930684
#> [3,] -0.17713234  0.7999603  -7.862937
#> [4,] -0.02946255  0.2816377  -5.191110
#> [5,]  0.06848906  0.3584364  -5.859788
#> [6,] -0.17410827  0.7329111  -2.059854
head(mbe_fit$weight_data[, c("clinical", "surrogate1", "surrogate2", "norm_w")])
#>       clinical surrogate1 surrogate2       norm_w
#>          <num>      <num>      <num>        <num>
#> 1: -0.10720750  0.7560167  -7.127875 0.0002672069
#> 2:  0.01877278  0.6573292  -7.300413 0.0002656093
#> 3: -0.11400997  0.8274257  -4.235286 0.0002099456
#> 4: -0.21420967  1.0618332  -4.698719 0.0002687302
#> 5: -0.04365231  0.1859734  -3.115996 0.0002484868
#> 6: -0.14840752  0.6815660  -3.572072 0.0002695261
```
