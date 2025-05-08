Gibbs sampling with block Metropolis
================
Hyejung Lee <hyejung.lee@utah.edu>
Wed May 07, 2025 18:27:31

- [Overview](#overview)
- [Setup](#setup)
- [Data](#data)

<!-- Now write your README content in Markdown, plus R chunks if you like: -->

# Overview

This directory contains `gibbs_sampler_v2.R` file, which is a function
that performs Gibbs sampling with block Metropolis for acute and chronic
GRF slopes. Use of Metropolis for regression coefficients really helps
to explore full prior supports because the regression coefficients are
highly correlated to each other.

  

# Setup

Let:

- $\hat{\theta}_i$ = (observed) treatment effect on the clinical
  endpoint from $i^{th}$ study
- $\hat{\gamma}_{1i}$ = (observed) treatment effect on the chronic slope
  from $i^{th}$ study
- $\hat{\gamma}_{2i}$ = (observed) treatment effect on the acute slope
  from $i^{th}$ study

First stage model

$$
\begin{bmatrix} 
\hat{\theta}_{i} \\\ 
\hat{\gamma}_{1i} \\\ 
\hat{\gamma}_{2i} 
\end{bmatrix} 
\bigg| 
\begin{bmatrix} 
\theta_{i} \\\ 
\gamma_{1i} \\\ 
\gamma_{2i}
\end{bmatrix} 
\sim 
N\left(
\begin{bmatrix} 
\theta_{i} \\\ 
\gamma_{1i} \\\ 
\gamma_{2i}
\end{bmatrix}, 
\begin{bmatrix}
\sigma^2_{i} & r_{\theta, 1i} \sigma_{i} \delta_{1i} & r_{\theta, 2i} \sigma_{i} \delta_{2i} \\\
r_{\theta, 1i} \sigma_{i} \delta_{1i} & \delta^2_{1i} & r_{1i, 2i} \delta_{1i} \delta_{2i} \\\
r_{\theta, 2i} \sigma_{i} \delta_{2i} & r_{1i, 2i} \delta_{1i} \delta_{2i} & \delta^2_{2i} \end{bmatrix} 
\right)$$

The second stage model:

$$
\begin{bmatrix} \theta_{i} \\ \gamma_{1i} \\ \gamma_{2i}\end{bmatrix}, \sim N\left(\begin{bmatrix} \mu_{\theta} \\ \mu_{\gamma 1} \\ \mu_{\gamma 2}\end{bmatrix}, \begin{bmatrix}
\sigma^2_{\theta} & R_{\theta, \gamma 1} \sigma_{\theta} \sigma_{\gamma 1} & R_{\theta, \gamma 2} \sigma_{i} \sigma_{\gamma 1} \\
R_{\theta, \gamma 1} \sigma_{\theta} \sigma_{\gamma 1} & \sigma^2_{\gamma 1} & R_{\gamma 1, \gamma 2} \sigma_{\gamma 1} \sigma_{\gamma 2}  \\
R_{\theta, \gamma 2} \sigma_{i} \sigma_{\gamma 2} & R_{\gamma 1, \gamma 2} \sigma_{\gamma 1} \sigma_{\gamma 2} & \sigma^2_{\gamma 2} 
\end{bmatrix}\right), i=1,...,m 
$$ {#eq:EQ2}

The second stage model can also be expressed in terms of linear
expression as below:

$E(\theta_i | \gamma_{1i}, \gamma_{2i}) = \alpha_\theta + \beta_1 \gamma_{1i} + \beta_2 \gamma_{2i}$

$Var(\theta_i | \gamma_{1i}, \gamma_{2i}) = \lambda^2_{\theta}$

$E(\gamma_{1i} | \gamma_{i2}) = \alpha_{\gamma_1} + \omega \gamma_{2i}$

$Var(\gamma_{1i} | \gamma_{2i}) = \lambda^2_{\gamma 1}$

$E(\gamma_{2i} ) = \mu_{\gamma2}$

$Var(\gamma_{2i} ) = \sigma^2_{\gamma2}$

  

The hyperparameters are called the following in the code:

- `alphaCEonSur1Sur2` $= \alpha_\theta$
- `b1CEonSur1Sur2` $= \beta_1$
- `b2CEonSur1Sur2` $= \beta_2$
- `SigSqCEonSur1Sur2` $= \lambda^2_{\theta}$
- `alphaSur1onSur2`$= \alpha_{\gamma_1}$
- `bSur1onSur2` $= \omega$
- `SigSqSur1onSur2` $=\lambda^2_{\gamma 1}$
- `muSur2` $=\mu_{\gamma2}$
- `sigSqSur2` $=\sigma^2_{\gamma2}$

  

# Data

First, you need to put your observed data into a list with the following
names:

- Sur1Est = estimated treatment effect on the chronic slope
- Sur1SE = standard error of treatment effect on the chronic slope
- Sur2Est = estimated treatment effect on the acute slope
- Sur2SE = standard error of treatment effect on the acute slope
- ClnEst = estimated treatment effect on the clinical endpoint
- ClnSE = standard error of treatment effect on the clinical endpoint
- R1Clin = correlation between treatment effect on clinical endpoint and
  chronic slope
- R2Clin = correlation between treatment effect on clinical endpoint and
  acute slope
- R1R2 = correlation between treatment effect on acute and chronic
  slopes

``` r
#Example data.
#suppose myData is the data frame that contains observed data of interest.
#create a list "dat_list" using this data frame with correct name:

dat_list<-data.table(
  Sur1Est = myData$beta31_beta32_est,
  Sur1SE  = myData$beta31_beta32_se,

  Sur2Est = myData$zbeta11_beta12_est,
  Sur2SE  = myData$zbeta11_beta12_se,

  ClnEst  = myData$loghr1_est,
  ClnSE   = myData$loghr1_se,

  R1Clin = myData$beta3_ce_corr,
  R2Clin = myData$beta1_ce_corr,
  R1R2   = myData$dgfr_chr_corr
)
```

  

You also need to give some initial values. For me, I decided on 4
initial values:

1.  Values same as Tom’s SAS initial values
2.  One standard deviation **smaller** than the mean of the prior
    distribution for the regression coefficients, and 25th percentile of
    the prior distribution for the variance
3.  Median of the prior distribution
4.  One standard deviation **bigger** than the mean of the prior
    distribution for the regression coefficients, and 75th percentile of
    the prior distribution for the variance

``` r
init_list<-list(
  #this is same as Tom's SAS code
  list(
    "b2CEonSur1Sur2"=0,
    "SigSqCEonSur1Sur2"=0.000010,
    "alphaCEonSur1Sur2"=0,
    "b1CEonSur1Sur2"=0,
    
    "SigSqSur1onSur2"=1,
    "alphaSur1onSur2"=0,
    "bSur1onSur2"=0,
    
    "muSur2"=-.4645,
    "sigSqSur2"=2.1993
  ),
  
  
  #this is at the lower value
  list(
    "alphaCEonSur1Sur2"=-100, #One standard deviation below
    "b2CEonSur1Sur2"=-100, #One standard deviation below
    "b1CEonSur1Sur2"=-100, #One standard deviation below
    "SigSqCEonSur1Sur2"=qinvgamma(0.25, shape=0.261,rate=0.000408), #this is the cut off value for low heterogeneity for SD
    
    "alphaSur1onSur2"=-100, #One standard deviation below
    "bSur1onSur2"=-100, #One standard deviation below
    "SigSqSur1onSur2"=qinvgamma(0.25, shape=0.261,rate=0.005),#this is the cut off value for low heterogeneity for SD
    
    "muSur2"=-100, #One standard deviation below
    "sigSqSur2"=qinvgamma(0.25, shape=0.261,rate=0.000408)#this is the cut off value for low heterogeneity
  ),
  
  
  #this is at the middle value
  list(
    "alphaCEonSur1Sur2"=0, #mean
    "b2CEonSur1Sur2"=0, #mean
    "b1CEonSur1Sur2"=0, #mean
    "SigSqCEonSur1Sur2"=qinvgamma(0.5, shape=0.261,rate=0.000408), #median
    
    "alphaSur1onSur2"=0,
    "bSur1onSur2"=0,
    "SigSqSur1onSur2"=qinvgamma(0.5, shape=0.261,rate=0.005),#this is the cut off value for low heterogeneity for SD
    
    "muSur2"=0, #One standard deviation below
    "sigSqSur2"=qinvgamma(0.5, shape=0.261,rate=0.000408)#this is the cut off value for low heterogeneity
  ),
  
  
  #this is at the higher value
  list(
    "alphaCEonSur1Sur2"= 100, #One standard deviation above
    "b2CEonSur1Sur2"= 100, #One standard deviation above
    "b1CEonSur1Sur2"= 100, #One standard deviation above
    "SigSqCEonSur1Sur2"=qinvgamma(0.75, shape=0.261,rate=0.000408), #this is the cut off value for low heterogeneity for SD
    
    "alphaSur1onSur2"= 100, #One standard deviation above
    "bSur1onSur2"= 100, #One standard deviation above
    "SigSqSur1onSur2"=qinvgamma(0.75, shape=0.261,rate=0.005),#this is the cut off value for medium heterogeneity for SD
    
    "muSur2"= 100, #One standard deviation above
    "sigSqSur2"=qinvgamma(0.75, shape=0.261,rate=0.000408) 
  )
  
)
```

We are now ready to run!

``` r
source("gibbs_sampler_v2.R")
set.seed(30)
tic() #If you want to keep track of time
gibbs_out2<-gibbs_sampler_v2(dat_list, #my data
                             inits = init_list, #list of initial values. Number of list is the number of chains
                             burn.out = 0.5, #proportion of n_iter as warm up
                             n_iter = 70000, #number of iterations
                             parallel = T, #parallel run for different chainas. 
                             ncores = 5) #number of cores to use for parallel. If there are lesser number of chains than the specified ncores, let's say 4 chains, then only 4 cores will be used.
toc()
#104.976 sec elapsed
```
