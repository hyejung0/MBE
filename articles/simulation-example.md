# Simulation Example

``` r

library(MBE)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
```

We use the simulated data provided in the package to demonstrate an MBE
with one clinical endpoint and exactly two surrogate endpoints. In this
CKD example, the surrogates are chronic and acute eGFR slopes. The
simulation scenarios are summarized in the documentation for
`interim_sim_dat`.

## Data

The simulated data is stored in the `trial_sim_dat` dataset in the
package. The dataset provides estimated treatment effects, standard
errors, and correlations on the clinical endpoint, chronic slope, and,
acute slope in 66 trials as well as the corresponding true treatment
effects. Below, we show the first few lines of the dataset followed by
column descriptions.:

``` r

data("trial_sim_dat")
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

| Column Name | Description |
|---:|:---|
| trial id | Integer. Index for trial. |
| CE_est | Numeric. Observed estimated treatment effect on the clinical endpoint. |
| Sur1_est | Numeric. Observed estimated treatment effect on the first surrogate endpoint. |
| Sur2_est | Numeric. Observed estimated treatment effect on the second surrogate endpoint. |
| CE_se | Numeric. Standard error of the observed estimated treatment effect on the clinical endpoint. |
| Sur1_se | Numeric. Standard error of the observed estimated treatment effect on the first surrogate endpoint. |
| Sur2_se | Numeric. Standard error of the observed estimated treatment effect on the second surrogate endpoint. |
| Cor_CE_Sur1 | Numeric. Estimated correlation between clinical endpoint and 1st surrogate endpoint. |
| Cor_CE_Sur2 | Numeric. Estimated correlation between clinical endpoint and 2nd surrogate endpoint. |
| Cor_Sur1_Sur2 | Numeric. Estimated correlation between 1st surrogate endpoint and 2nd surrogate endpoint. |
| theta_CE | Numeric. True treatment effect on the clinical endpoint. |
| theta_Sur1 | Numeric. True treatment effect on the 1st surrogate endpoint. |
| theta_Sur2 | Numeric. True treatment effect on the 2nd surrogate endpoint. |

For our demonstration, we will just take out the first trial as if it is
a ‘new trial’ and fit the historical modeling on the remaining 65 trial.
We will then estimate MBE on the trial that was left out.

## Historical Model Fitting

Here, we demonstrate how to fit the historical model using the stan
script `random_beta0_invGamma.stan` which allows $`\beta_0`$ to be
estimated from the data and uses inverse gamma prior for the variance
parameters.

First, we reshape the data into a list of observed treatment effects and
covariance matrices for each study. The observed treatment effects are
stored in a list of vectors, where each vector contains the observed
treatment effects on the clinical endpoint, chronic slope, and acute
slope for a given study. The observed covariance matrices are stored in
a list of matrices, where each matrix contains the variances and
covariances of the observed treatment effects for a given study.

``` r


library(cmdstanr)


# Generate data for Stan model fitting. The data is a list of the following elements:

# observed treatment effect per study
obs_mean<-lapply(1:nrow(trial_sim_dat),function(j){
  unlist(trial_sim_dat[j,.(CE_est,      ## observed Clinical endpoint
                    Sur1_est,     ## observed chronic slope
                    Sur2_est      ## observed acute slope
  )])
})


# observed covariance matrix per study
obs_var <- lapply(1:nrow(trial_sim_dat), function(j) {
  # vector of SEs
  sample_dat<-trial_sim_dat[j,.(CE_se, #SE for Clinical endpoint
                         Sur1_se, #SE for chronic slope
                         Sur2_se, #SE for acute slope
                         Cor_CE_Sur1, #correlation of chronic slope & clinical endpoint
                         Cor_CE_Sur2, #correlation of acute slope & clinical endpoint
                         Cor_Sur1_Sur2  #correlation of acute & chronic slope
  )]


  sample_dat[
    ,
    matrix(
      c(CE_se^2, Cor_CE_Sur1*CE_se*Sur1_se,  Cor_CE_Sur2*CE_se*Sur2_se,
        Cor_CE_Sur1*CE_se*Sur1_se, Sur1_se^2, Sur1_se*Sur2_se*Cor_Sur1_Sur2,
        Cor_CE_Sur2*CE_se*Sur2_se, Sur1_se*Sur2_se*Cor_Sur1_Sur2,Sur2_se^2),
      nrow=3,ncol = 3, byrow=T)
  ]
})


#collect them together as a more comprehensive list
#Exclude the first RCT from the list
data_sub = list(
  N=nrow(trial_sim_dat) -1,
  obs_mean = obs_mean[-1],
  obs_var = obs_var[-1]

)
```

Now that we have the data in the correct format, we can fit the
historical model using the `cmdstanr` package. We will use the
`random_beta0_invGamma.stan` script to fit the model, which allows for
$`\beta_0`$ to be estimated from the data and uses an inverse gamma
prior for the variance parameters.

To fit the model, we use our MCMC function on the sample data. *Note:
This step may take several minutes to run on your machine.*

``` r


# Read in the stan code: historical model with inverse gamma prior for variance parameters and model $\beta_0$ as a random parameter
mod_path <- system.file("stan", "random_beta0_invGamma.stan", package = "MBE")
cat(readLines(mod_path), sep = "\n")
#> // One clinical endpoint + exactly two surrogate endpoints.
#> // Estimate beta0 and use inverse-gamma priors for variances.
#> 
#> data {
#>   int<lower=1> N;                        // number of studies
#>   array[N] vector[3] obs_mean;           // fixed order: CE, surrogate 1, surrogate 2
#>   array[N] matrix[3, 3] obs_var;         // covariance for those three effects
#> }
#> parameters{
#> 
#>   //Population mean
#>   real muSur2;              // mean for acute slope
#> 
#>   //meta-regression for modeling chronic slope ~ acute slope
#>   real alphaSur1onSur2;     // intercept
#>   real bSur1onSur2;         // slope
#> 
#> //meta-regression for modeling clinical endpoints ~ chronic slope + acute slope
#>   real alphaCEonSur1Sur2;   // intercept (= beta0)
#>   real b1CEonSur1Sur2;      // slope for chronic slope
#>   real b2CEonSur1Sur2;      // slope for acute slope
#> 
#>   // Variance components
#>   real<lower=0> sigSqSur2;          //variance for acute slope
#>   real<lower=0> SigSqSur1onSur2;    //conditional variance for modeling chronic slope
#>   real<lower=0> SigSqCEonSur1Sur2;  //conditional variance for modeling Clinical endpoint
#> 
#>   // standard normal distribution to construct psi later
#>   matrix[3, N] z;
#> }
#> 
#> transformed parameters {
#> 
#>   // Construct 'psi' here using the non-centered parameterization
#>   array[N] vector[3] psi;  // Latent values per study (now an array of vectors)
#> 
#> 
#>   // Construct psi using the reparameterization formula: psi = mu + sqrt(variance) * z
#>   // We loop through each study (column of z)
#>   for (i in 1:N) {
#> 
#>     psi[i][3]=muSur2+sqrt(sigSqSur2)*z[3,i];
#>     psi[i][2]=alphaSur1onSur2 + bSur1onSur2*psi[i][3]+sqrt(SigSqSur1onSur2)*z[2,i];
#>     psi[i][1]=alphaCEonSur1Sur2 + b1CEonSur1Sur2*  psi[i][2]   + b2CEonSur1Sur2* psi[i][3]+sqrt(SigSqCEonSur1Sur2)*z[1,i];
#>   }
#> }
#> 
#> 
#> 
#> 
#> model {
#>   // Priors
#>   muSur2 ~ normal(0, 100);
#>   sigSqSur2~ inv_gamma(0.261, 0.005);
#> 
#>   alphaSur1onSur2 ~ normal(0, 100);
#>   bSur1onSur2 ~ normal(0, 100);
#>   SigSqSur1onSur2~ inv_gamma(0.261, 0.005);
#> 
#>   alphaCEonSur1Sur2 ~ normal(0, 100);
#>   b1CEonSur1Sur2 ~ normal(0, 100);
#>   b2CEonSur1Sur2 ~ normal(0, 100);
#>   SigSqCEonSur1Sur2 ~ inv_gamma(0.261, 0.000408);
#> 
#> 
#>   to_vector(z) ~ std_normal();
#> 
#>   // The Stage 1 model (psi ~ multi_normal(mu, Sigma)) is now implicitly
#>   // defined by the construction of psi in the transformed parameters block.
#> 
#>   // Stage 2 model (Likelihood) - this remains the same
#>   for (i in 1:N) {
#>     obs_mean[i] ~ multi_normal(psi[i], obs_var[i]);
#>   }
#> 
#> }
#> 
#> generated quantities {
#>   vector[N] log_lik;
#>   for (i in 1:N) {
#>     // The log-likelihood of the observed data. Used for leave-one-out cross-validation or other model comparison metrics.
#>     log_lik[i] = multi_normal_lpdf(obs_mean[i] | psi[i], obs_var[i]);
#>   }
#> 
#> 
#>   //Derive marginal variance for Sur1 and CE
#>   real sigSqSur1 = SigSqSur1onSur2 + square(bSur1onSur2)*sigSqSur2;
#>   real sigSqClin = SigSqCEonSur1Sur2 +
#>                     square(b1CEonSur1Sur2)*sigSqSur1 +
#>                     square(b2CEonSur1Sur2)*sigSqSur2 +
#>                     2.0 * b1CEonSur1Sur2 * b2CEonSur1Sur2 * bSur1onSur2 * sigSqSur2;
#> 
#>   //R^2 for regression on CE~sur1 + sur2
#>   real R2CEonSur1Sur2 = 1.0 - (SigSqCEonSur1Sur2/sigSqClin);
#> 
#>   //RMSE for regression on CE~sur1 + sur2
#>   real RMSECEonSur1Sur2 = sqrt(SigSqCEonSur1Sur2);
#> 
#> }
```

``` r



# Run the MCMC sampling

ncores<-4 #Number of cores to run chains parallely
options(mc.cores = ncores) #necessary for setting parallel
num.iter=40000 #number of iterations for each chain. You can change this number as you like. We recommend to have at least 40,000 iterations for the historical model fitting.
num.warmup=num.iter/2 #number of warmup iterations. The default is half of the total number of iterations, but you can change this number as you like. We recommend to have at least 20,000 warmup iterations for the historical model fitting.
num.chains<-ncores
num.thin<-1 #thinning number. The default is 1, which means no thinning. You can change this number as you like. We recommend to have no thinning for the historical model fitting.


set.seed(1) #Change seed number as you like

fit <- mod$sample(
  data = data_sub,
  chains = num.chains,
  parallel_chains = num.chains,
  iter_warmup = num.warmup,
  iter_sampling = num.warmup,
  adapt_delta = 0.99,
  max_treedepth=25,
  output_dir = ".",  #current working directory is where the STAN output will be saved.
  output_basename="historical_fit" #the STAN files are saved with "historical_fit" as a base name
)
```

    #> All 4 chains finished successfully.
    #> Mean chain execution time: 118.0 seconds.
    #> Total execution time: 144.6 seconds.
    #> 149.839 sec elapsed

Inspect the posterior distribution of all parameters in the historical
model.:

``` r

all_pars<-c(
  "alphaCEonSur1Sur2" ,
  "b1CEonSur1Sur2",
  "b2CEonSur1Sur2",
  "SigSqCEonSur1Sur2",
  "alphaSur1onSur2",
  "bSur1onSur2",
  "SigSqSur1onSur2",
  "muSur2",
  "sigSqSur2"
)
#first, we can draw the historical posterior samples if we'd like
historical_posterior<-fit$draws(variables = all_pars, format = "draws_df")


# Generate the summary locally
summary_df <- fit$summary(variables = all_pars)
```

``` r

# Load and print the tiny summary table
print(summary_df)
#> # A tibble: 9 × 10
#>   variable         mean   median      sd     mad       q5     q95  rhat ess_bulk
#>   <chr>           <dbl>    <dbl>   <dbl>   <dbl>    <dbl>   <dbl> <dbl>    <dbl>
#> 1 alphaCEonSu… -0.0671  -0.0670  0.0482  0.0471  -1.46e-1  0.0119  1.00   35932.
#> 2 b1CEonSur1S… -0.295   -0.296   0.0639  0.0624  -3.98e-1 -0.189   1.00   34155.
#> 3 b2CEonSur1S… -0.0353  -0.0353  0.00495 0.00489 -4.33e-2 -0.0271  1.00   46408.
#> 4 SigSqCEonSu…  0.00404  0.00269 0.00411 0.00283  3.18e-4  0.0123  1.00   26435.
#> 5 alphaSur1on…  0.643    0.643   0.0880  0.0868   5.00e-1  0.788   1.00   16787.
#> 6 bSur1onSur2   0.0232   0.0233  0.0159  0.0156  -3.13e-3  0.0491  1.00   15326.
#> 7 SigSqSur1on…  0.254    0.242   0.0794  0.0733   1.46e-1  0.401   1.00   22123.
#> 8 muSur2       -2.49    -2.49    0.772   0.766   -3.77e+0 -1.24    1.00    9379.
#> 9 sigSqSur2    30.0     29.2     6.47    6.15     2.09e+1 41.8     1.00   14417.
#> # ℹ 1 more variable: ess_tail <dbl>
```

## Estimating Posterior Distribution of MBE

Now that we have historical posteriors, we use the `MBE` function to
estimate the posterior distribution of MBE for a new trial. We will use
the first trial in the `trial_sim_dat` dataset as our new trial and use
the remaining trials to fit the historical model.

``` r

# Put the new trial into the appropriate format

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


#Estimate MBE
MBE_distribution<-MBE(
mcmc_dat = historical_posterior,
sample_dat=one_sim_dat,
diffuse_se = 100,
diffuse = TRUE,
intercept0 = TRUE
)
```

In the above [`MBE()`](https://hyejung0.github.io/MBE/reference/MBE.md)
function, we set `intercept0=TRUE`, which centers the distribution of
$`\beta_0`$ at zero (the `alphaCEonSur1Sur2` column in
`historical_posterior`). We use a diffuse prior for the two surrogate
effects (`diffuse=TRUE` and `diffuse_se=100`) to avoid strong prior
information about those effects. The diffuse prior is bivariate normal
with mean zero and marginal standard deviation 100. The
[`MBE()`](https://hyejung0.github.io/MBE/reference/MBE.md) function
returns the posterior distribution stored in `MBE_distribution`.

The posterior distribution is trivariate normal with mean and variance
as below:

``` r

# print out mean
MBE_distribution$post_mean
#>                  [,1]
#> clinical   -0.1013294
#> surrogate1  0.6837417
#> surrogate2 -5.4760921

# print out variance
MBE_distribution$post_var
#>               clinical   surrogate1   surrogate2
#> clinical    0.01196417 -0.019592518 -0.091338433
#> surrogate1 -0.01959252  0.060139164 -0.009146439
#> surrogate2 -0.09133843 -0.009146439  3.981621682
```

Thus, the posterior distribution of MBE on the treatment effect on the
clinical endpoint is normal with mean -0.101 and variance 0.012.
Similarly, the posterior distribution of MBE on the treatment effect on
the chronic slope is \$N(\$0.684,0.060 $`)`$. The posterior distribution
of MBE on the treatment effect on the acute slope $`N(`$-5.476, 3.982
$`)`$.

The posterior distribution samples can also be obtained with
`MBE_distribution$post_psi0`, which has the same number of rows as the
number of MCMC samples in the historical model fitting. Each row is a
sample from the posterior distribution of MBE, and each column
corresponds to the treatment effect on the clinical endpoint, chronic
slope, and acute slope, respectively.

``` r

head(MBE_distribution$post_psi0)
#>         clinical surrogate1 surrogate2
#> [1,] -0.01742418  0.4358344  -7.647448
#> [2,] -0.05070114  0.3625249  -5.290546
#> [3,] -0.21549799  0.9284195  -4.653444
#> [4,] -0.22270912  1.0586468  -4.846612
#> [5,] -0.02885444  0.6432777  -8.217527
#> [6,] -0.22661395  0.9851503  -4.977536
```

The `MBE_distribution$weight_data` is a data.table which contains
unweighted posterior draws and the associated weights for those draws.
The `MBE_distribution$weight_data` has the following columns.:

``` r

head(MBE_distribution$weight_data)
#>       clinical surrogate1 surrogate2       norm_w         w
#>          <num>      <num>      <num>        <num>     <num>
#> 1: -0.08741164  0.8217200 -10.292497 0.0002672069 0.8190106
#> 2: -0.21601529  0.8178430  -3.292377 0.0002656093 0.8141140
#> 3: -0.20075066  0.7338019  -5.747120 0.0002099456 0.6435002
#> 4: -0.07627662  0.5728434  -4.210812 0.0002687302 0.8236798
#> 5:  0.12683662  0.1595884  -4.466493 0.0002484868 0.7616320
#> 6: -0.18528890  0.7034575  -5.451308 0.0002695261 0.8261194
```

Please use the `norm_w` (normalized weight) to obtain the proper
posterior distribution. The `w` is unnormalized weight.
