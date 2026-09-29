# Interim Analysis Example

``` r

library(MBE)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
```

We use the simulated data provided in the package to demonstrate the
input structure for an interim MBE analysis. The packaged data contain
one clinical endpoint and exactly two surrogate endpoints: chronic
(`Sur1`) and acute (`Sur2`) eGFR slopes. No third surrogate is used. The
simulation scenarios are summarized in the documentation for
`interim_sim_dat`.

``` r

data("interim_sim_dat")
names(interim_sim_dat)
#>  [1] "case"           "analysis.month" "R1Clin"         "R2Clin"        
#>  [5] "R12"            "ClnEst"         "ClnSE"          "Sur2Est"       
#>  [9] "Sur2SE"         "Sur1Est"        "Sur1SE"
```
