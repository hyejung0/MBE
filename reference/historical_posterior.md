# Historical Posterior for Two Surrogate Endpoints

Posterior draws of all model parameters from the MCMC model with chronic
and acute eGFR slopes as the two surrogates for the clinical endpoint.
This data set is used in the examples of the package to demonstrate how
to use the `vector_to_matrix` and `MBE` functions.

## Usage

``` r
historical_posterior
```

## Format

A data.table with 4,000 rows and 9 variables:

- alphaCEonSur1Sur2:

  Intercept for clinical effect conditional on both surrogate effects.

- b1CEonSur1Sur2:

  Slope for surrogate 1 in the clinical-effect regression.

- b2CEonSur1Sur2:

  Slope for surrogate 2 in the clinical-effect regression.

- SigSqCEonSur1Sur2:

  Conditional variance of the clinical effect.

- alphaSur1onSur2:

  Intercept for surrogate 1 conditional on surrogate 2.

- bSur1onSur2:

  Slope for surrogate 2 in the surrogate-1 regression.

- SigSqSur1onSur2:

  Conditional variance of surrogate 1.

- muSur2:

  Population mean for surrogate 2.

- sigSqSur2:

  Population variance for surrogate 2.
