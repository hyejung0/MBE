
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Multi-component Bayesian Endpoint (MBE)

<!-- badges: start -->

<!-- badges: end -->

This repository contains all relevant R codes that was used in the
Multi-Component Bayesian Endpoints for Randomized Clinical Trials paper
produced by Lee et. al., (2026).

The MBE is an endpoint that combines clinical endpoint with at least one
surrogate endpoint to create a more sensitive and efficient endpoint for
clinical trials. It is developed using a Bayesian framework, which
allows for the incorporation of historical randomized clinical trials
(RCTs) and the estimation of uncertainty in estimating the treatment
effect on the clinical endpoint using the treatment effect on surrogate
endpoint(s).

Thus, it has two distinct coding components:

- Step 1. Fit Markov chain Monte Carlo (MCMC) on a trial-level model to
  historical RCTs to produce historical posterior samples of the model
  parameters

- Step2. Construct a prior distribution of MBE using the historical
  posterior samples estimated in step 1, and then update the MBE’s prior
  based on the estimated treatment effects in a “new RCT” to obtain the
  MBE’s posterior distribution for the true treatment effect on the CE
  in that new RCT.

For both the historical and new RCTs, the model is a hierarchical model
that relates the treatment effect on the clinical endpoint to the
treatment effect on surrogate endpoint(s) across trials. The model
parameters include the intercept, slope(s), and variance parameters that
capture the relationship between the treatment effects on the surrogate
and clinical endpoints across trials. Please refer to Section 2 of the
manuscript for details of the model specification.

Our work demonstrated these steps using chronic kidney disease trials
where two surrogate endpoints, the treatment effect on acute and chronic
slopes of glomerular filtration rate (GFR), were used to estimate the
MBE for the treatment effect on the clinical endpoint.

## Data Structure

We have 66 historical RCTs that have treatment effect estimates on the
clinical endpoint, the acute slope, and the chronic slope. Data use
agreements prohibit CKD-EPI CT from sharing data with parties external
to the agreement. Thus, we explain the data structure and how to use the
data to fit the historical model, but we cannot share the data itself.

As the analysis is at trial-level, the unit of observation is per trial.
Thus, the dataset has 66 rows. The columns are as follows:

| Column Name | Description |
|----|----|
| ClnEst | observed treatment effect on Clinical endpoint (CE) |
| ClnSE | standard error of observed treatment effect on CE |
| Sur1Est | observed treatment effect on chronic slope (Surrogate 1) |
| Sur1SE | standard error of observed treatment effect on chronic slope |
| Sur2Est | observed treatment effect on acute slope (Surrogate 2) |
| Sur2SE | standard error of observed treatment effect on acute slope |
| R1Clin | correlation between observed treatment effect on CE and chronic slope |
| R2Clin | correlation between observed treatment effect on CE and acute slope |
| R1R2 | correlation between observed treatment effect on chronic slope and acute slope |

These column names are used in all R codes provided in this Github. So
if you want to use the provided R script to your own data, please make
sure to have the same column names and structure as described above.

## Step 1 Historical Model Fitting

The manuscript describes two different models to fit the historical
data:

1.  Fix $\beta_0$ to be zero
2.  Allow $\beta_0$ to be estimated from the data

where $\beta_0$ is the intercept parameter in the regression model that
relates the true treatment effect on the clinical endpoint to the true
treatment effect on all surrogate endpoints.

Besides, we can also choose to use either inverse gamma or half-normal
priors for the variance parameters (or standard deviation) in the model.
Thus, totally, there are four different ways to fit the historical
model.

The four stan code that correspond to the four different ways to fit the
historical model are stored in the `code/historical` folder. :

- `fix_beta0_invGamma.stan`: fit the historical model with inverse gamma
  prior for variance parameters and $\beta_0$ fixed to be zero
- `random_beta0_invGamma.stan`: fit the historical model with inverse
  gamma prior for variance parameters and model $\beta_0$ as a random
  parameter
- `fix_beta0_halfNormal.stan`: fit the historical model with half-normal
  prior for standard deviation parameters and $\beta_0$ fixed to be zero
- `random_beta0_halfNormal.stan`: fit the historical model with
  half-normal prior for standard deviation parameters and model
  $\beta_0$ as a random parameter

In the same folder, you can find `historical_model_fit.R` file, which
sources one of the four stan files to run MCMC on the historical data.
You can choose which stan file to source based on the model
specification you want to fit.

## Step 2 MBE Construction

Using the parameter distributions estimated with the historical data in
Step 1, we now estimate the distribution of MBE. All relevant R code for
this step can be found in the `code/MBE` folder.

The function `theta0_given_hat_psi0` estimates the posterior
distribution of $\theta_0 \mid \hat{\psi}_0$ (saved as
`theta0_given_hat_psi0.R` file). This function fits Models A through D
in the manuscript. `MBE_given_all_dat.R` file demonstrates use of this
function to generate the posterior distributions of
$\theta_0 \mid \hat{\psi}_0$ under all four different model
specifications.

with full historical prior carry over is shown in
`MBE_full_historical_posterior.R` file. This is what the manuscript
calls ‘Model A.’

The R code for Step 2 is stored in the `code/MBE_estimation.R` file. In
this code, we first read in the MCMC samples of the model parameters
estimated in Step 1, and then we use these samples to construct the MBE
for a “new RCT” that was left out from the historical data. The “new
RCT” is represented by the `sample_dat` list in the code, which contains
the observed treatment effect estimates on the clinical endpoint and
surrogate endpoints, as well as their standard errors and correlations.

## Historical Data Model Fit Check

In Appendix, you can find a section for Historical Model fit, where I
used a various metric to check the model fit for the historical data.
The R code for this is stored in the `code/historical_model_assessment`
folder. Suppose that as in our paper, we have 66 historical RCTs.

### Leave-One-Out Cross-Validated Tolerance Interval (LOO-CV TI) coverage

The `PPD.R` file demonstrates how to construct a posterior predictive
distribution for observed treatment effect on the clinical endpoint. One
fits historical model on 65 RCTs, and then uses the posterior samples of
the model parameters to generate a posterior predictive distribution
(PPD) for the observed treatment effect on the clinical endpoint for the
left out RCT. The PPD is estimated by following procedure. Suppose we
have $B$ posterior samples of the model parameters, and we have $N$ RCTs
in the historical data. Then for each left out RCT, we do the following:

1.  Estimate the posterior distribution of the surrogate endpoint(s) by
    updating the historical model with the observed treatment effect on
    the surrogate endpoint(s) for the left out RCT. This should be done
    for each of the $B$ posterior samples of the model parameters and we
    will have $B$ posterior samples of the surrogate endpoint(s) for the
    left out RCT.
2.  For each $b$ posterior sample of the surrogate endpoint(s), we
    generate a posterior predictive distribution for the observed
    treatment effect on the clinical endpoint for the left out RCT by
    drawing a sample from each conditional normal distribution:

$$N(\beta_0^{(b)} + \beta_1^{(b)} \cdot \gamma_{01}^{(b)}+ \beta_2^{(b)} \cdot \gamma_{02}^{(b)}, \lambda_{\theta}^{2 \ (b)} + \sigma^2_{\hat{\theta}_0}),$$

where $\beta_0^{(b)}$, $\beta_1^{(b)}$, $\beta_2^{(b)}$, and
$\lambda_{\theta}^{2 \ (b)}$ are $b^{th}$ MCMC sample of the the
historical posterior distribution, $\gamma_{01}^{(b)}$ and
$\gamma_{02}^{(b)}$ are the $b^{th}$ posterior samples of the surrogate
endpoint(s) for the left out RCT, and $\sigma^2_{\hat{\theta}_0}$ is the
observed variance of the treatment effect on the clinical endpoint for
the left out RCT.

The collection these draws form a posterior predictive distribution for
the observed treatment effect on the clinical endpoint for the left out
RCT. Then, we estimate where the observed treatment effect on the CE for
the left out RCT ($\hat{\theta}_0$) falls within this PPD. If we are
interested in constructing 95% tolerance interval, we can find whether
$\hat{\theta}_0$ sits within the 2.5th and 97.5th percentiles of this
PPD. We log the result as 1 if $\hat{\theta}_0$ falls within the 95%
tolerance interval, and 0 otherwise.

We repeat this process for all $N$ RCTs in the historical data, and
calculate the coverage of the 95% tolerance interval by taking the
average of the logged results (1 if $\hat{\theta}_0$ falls within the
95% tolerance interval, and 0 otherwise) across all $N$ RCTs. This
coverage should be close to 0.95 if the historical model fits the data
well.
