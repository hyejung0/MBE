# Simulated Data for Interim Analysis

This is a simulated data set for chronic kidney disease clinical trials.
The definitive clinical endpoint is the first of ESRD or a 57% decline
from baseline eGFR. In addition to the treatment effect estimated on the
definitive clinical endpoint, the treatment effect is estimated on
exactly two surrogate endpoints: chronic eGFR slope (`Sur1`) and acute
eGFR slope (`Sur2`).

The data contain eight scenarios. Cases 1 and 2 use a baseline eGFR
range of 30–60, cases 3 and 4 use 30–75, and cases 5–8 use 25–50
ml/min/1.73m2. Odd-numbered cases have no slope treatment effect;
even-numbered cases have an active slope treatment effect of 0.8. Cases
7 and 8 repeat the settings of cases 5 and 6, respectively.

## Usage

``` r
interim_sim_dat
```

## Format

A data frame with 528 rows and 11 variables:

- case:

  Integer. The scenario case number (1 through 8).

- analysis.month:

  Integer. Number of months since the study began for the specific case.

- R1Clin:

  Numeric. Estimated correlation between chronic slope and clinical
  endpoint. If NA, that is because the treatment effect on clinical
  endpoint (ClnEst and ClnSE) were not estimated due to small event
  size.

- R2Clin:

  Numeric. Estimated correlation between acute slope and clinical
  endpoint. If NA, that is because the treatment effect on clinical
  endpoint (ClnEst and ClnSE) were not estimated due to small event
  size.

- R12:

  Numeric. Estimated correlation between chronic slope and acute slope.

- ClnEst:

  Numeric. Estimated treatment effect on the clinical endpoint (log
  hazard ratio). If NA, that is because the treatment effect on clinical
  endpoint (ClnEst and ClnSE) were not estimated due to small event
  size.

- ClnSE:

  Numeric. Standard error of the estimated treatment effect on the
  clinical endpoint. If NA, that is because the clinical treatment
  effect was not estimated due to small event size.

- Sur2Est:

  Numeric. Estimated treatment effect on acute slope.

- Sur2SE:

  Numeric. Standard error of the estimated treatment effect on acute
  slope.

- Sur1Est:

  Numeric. Estimated treatment effect on chronic slope.

- Sur1SE:

  Numeric. Standard error of the estimated treatment effect on chronic
  slope.

## Source

Simulated by Jian Ying

## Author

Jian Ying <jian.ying@hsc.utah.edu>
