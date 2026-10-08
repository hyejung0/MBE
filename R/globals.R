#Here, we import all dependencies and set up global variables for the package.

#' @import data.table
NULL


#The Issue: R CMD check inspects code for variables that look like
#unquoted global variables (., ..these_cols, empirical_cdf, Sigma11, etc.)
#used inside data.table syntax in functions like vector_to_matrix().
utils::globalVariables(c(
  ".",
  "..these_cols",
  "empirical_cdf",
  "clinical",
  "surrogate1",
  "surrogate2",
  "norm_w",
  "w",
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
  "mu3",
  "det_Sigma22_33",
  "inv_Sigma22",
  "inv_Sigma23",
  "inv_Sigma33",
  "hat_psi02",
  "hat_psi03",
  "hat_Sigma_y0_22",
  "hat_Sigma_y0_23",
  "hat_Sigma_y0_33",
  "entry22",
  "entry33",
  "entry23",
  "inv_entry22",
  "inv_entry23",
  "inv_entry33",
  "post_Sigma22",
  "post_Sigma23",
  "post_Sigma33",
  "post_mu2",
  "post_mu3",
  "post_CE",
  "post_psi2",
  "post_psi3"
))
