
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MBE

<!-- badges: start -->

<!-- badges: end -->

This repository contains all relevant R functions that estimates
Multicomponent Bayesian Endpoints (MBE).

## Code for acute and chronic joint slopes

The `code/AcuChron` folder contains R scripts and data to accomplish two
objectives for developing the slope-based MBE based on the established
clinical endpoint (CE) and the acute and chronic GFR slopes. There are
two primary steps:

1.  Fit a trial-level model to historical RCTs to produce MCMC samples
    from the posterior distribution for the treatment effects on each of
    the 3 endpoints (CE and the acute and chronic GFR slopes)
2.  Construct a prior distribution linking the treatment effects on the
    3 endpoints from step 1, and then update the prior based on
    estimated treatment effects in a “new RCT” to obtain the postrior
    distribution for the true treatment effect on the CE in that new RCT

Below shows the explanation of each file/folder.

- `run_MCMC_repar.R`: Step 1 R script that runs MCMC on the historical
  trial-level data.
  - dependent files & folders:
    - `evt3data.rds`: Data from historical RCTs used to fit the trial
      level model for step 1 above. This data set contains the estimated
      treatment effects on the CE and slope-based endpoints across the
      66 EMA studies along with their associated standard errors and
      correlations
    - `./run_MCMC_output`: a folder that saves the STAN output csv file
      and the MCMC posterior samples from the trial level model
- `MBE_calc.R`: Step 2 R script that returns information concerning the
  posterior distribution of the true treatment effect on the CE in a new
  RCT
  - dependent files/folders:
    - `evt3data.rds`: Default inputs from the new RCT, which is set to
      be one of the previous 66 RCTs. This can be modified as desired to
      represent the new RCT of interest.
    - `./run_MCMC_output/ALL_STUDIES.rds`:posterior samples saved from
      the Step 1 rscript run_MCMC.R
    - `post_theta0_samples_v4.R`: R script that contains all functions
      required by `MBE_calc.R`

## Code for acute and chronic joint slopes using Maximum Likelihood Estimation

The `code/AcuChron` folder contains R script and file necessary to fit
maximum likelihood on hierarchical model for treatment effect on the CE
on treatment effect on acute and chronic slopes.

- `MLE.R`: fits MLE using all 66 studies
