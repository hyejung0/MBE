#Here, we import all dependencies and set up global variables for the package.

#' @import data.table
#' @import parallel
NULL


#The Issue: R CMD check inspects code for variables that look like
#unquoted global variables (., ..these_cols, empirical_cdf, Sigma11, etc.)
#used inside data.table syntax in functions like vector_to_matrix().
utils::globalVariables(c(
  ".",
  "..these_cols",
  "empirical_cdf",
  "bSur1onSur2",
  "SigSqSur1onSur2",
  "sigSqSur2",
  "muSur2",
  "alphaSur1onSur2",
  "alphaCEonSur1Sur2",
  "Sigma11",
  "Sigma12",
  "Sigma13",
  "Sigma22",
  "Sigma23",
  "Sigma33",
  "b1CEonSur1Sur2",
  "b2CEonSur1Sur2",
  "SigSqCEonSur1Sur2",
  "mu1",
  "mu2",
  "mu3"
))
