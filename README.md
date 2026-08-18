
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Multi-component Bayesian Endpoints (MBEs)

<!-- badges: start -->

[![R-CMD-check](https://github.com/hyejung0/MBE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hyejung0/MBE/actions/workflows/R-CMD-check.yaml)
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

Detailed explanation of mapping the manuscript to the code can be found
in this document: [Mapping the manuscript to the
code](articles/MBE-manuscript-to-code.html)

## Installation

You can install the development version of MBE from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("hyejung0/MBE")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MBE)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" alt="" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
