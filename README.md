
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MBE

<!-- badges: start -->

<!-- badges: end -->

This repository contains all relevant R functions that estimates
Multicomponent Bayesian Endpoints (MBE).

## Code for acute and chronic joint slopes

The `code/AcuChron` folder contains R scripts and data to accomplish two
objectives for developing the slope-based MBE based on the established
clinical endpoint (CE) and the acute and chronic GFR slopes.

The code is fixed to use inverse gamma priors for the variance
parameters when running MCMC on the historical data.

Please view the [README](./code/AcuChron/README.md) file in that folder
for details of each file and instructions on how to run the code.

## Code for acute and chronic joint slopes using Maximum Likelihood Estimation

The `code/AcuChron` folder contains R script and file necessary to fit
maximum likelihood on hierarchical model for treatment effect on the CE
on treatment effect on acute and chronic slopes.

- `MLE.R`: fits MLE using all 66 studies
