#' Historical Posterior for Two Surrogate Endpoints
#'
#' Posterior draws of all model parameters from the MCMC model with chronic
#' and acute eGFR slopes as the two surrogates for the clinical endpoint.
#' This data set is extracted from the `historical_fit` object and is provided
#' to avoid requiring users to refit the historical Stan model.
#' @format A data.table with 4,000 rows and 9 variables:
#' \describe{
#'   \item{alphaCEonSur1Sur2}{Intercept for clinical effect conditional on both surrogate effects.}
#'   \item{b1CEonSur1Sur2}{Slope for surrogate 1 in the clinical-effect regression.}
#'   \item{b2CEonSur1Sur2}{Slope for surrogate 2 in the clinical-effect regression.}
#'   \item{SigSqCEonSur1Sur2}{Conditional variance of the clinical effect.}
#'   \item{alphaSur1onSur2}{Intercept for surrogate 1 conditional on surrogate 2.}
#'   \item{bSur1onSur2}{Slope for surrogate 2 in the surrogate-1 regression.}
#'   \item{SigSqSur1onSur2}{Conditional variance of surrogate 1.}
#'   \item{muSur2}{Population mean for surrogate 2.}
#'   \item{sigSqSur2}{Population variance for surrogate 2.}
#' }
"historical_posterior"


#' Simulated Trial-Level Data for Two Surrogate Endpoints
#'
#' @description Simulated summary data for 66 randomized trials. Each trial has
#' one clinical endpoint and exactly two surrogate endpoints (chronic and acute
#' eGFR slopes), along with their standard errors, correlations, and simulation
#' truths.
#' @format A data frame with 66 rows (one row per trial) and 13 variables:
#' \describe{
#'   \item{trial_id}{Integer. Index for trial.}
#'   \item{CE_est}{Numeric. Observed estimated treatment effect on the clinical endpoint.}
#'   \item{Sur1_est}{Numeric. Observed estimated treatment effect on the first surrogate endpoint.}
#'   \item{Sur2_est}{Numeric. Observed estimated treatment effect on the second surrogate endpoint.}
#'   \item{CE_se}{Numeric. Standard error of the observed estimated treatment effect on the clinical endpoint.}
#'   \item{Sur1_se}{Numeric. Standard error of the observed estimated treatment effect on the first surrogate endpoint.}
#'   \item{Sur2_se}{Numeric. Standard error of the observed estimated treatment effect on the second surrogate endpoint.}
#'   \item{Cor_CE_Sur1}{Numeric. Estimated correlation between clinical endpoint and 1st surrogate endpoint.}
#'   \item{Cor_CE_Sur2}{Numeric. Estimated correlation between clinical endpoint and 2nd surrogate endpoint.}
#'   \item{Cor_Sur1_Sur2}{Numeric. Estimated correlation between 1st surrogate endpoint and 2nd surrogate endpoint.}
#'   \item{theta_CE}{Numeric. True treatment effect on the clinical endpoint.}
#'   \item{theta_Sur1}{Numeric. True treatment effect on the 1st surrogate endpoint.}
#'   \item{theta_Sur2}{Numeric. True treatment effect on the 2nd surrogate endpoint.}
#' }
#' @author Yizhen Xu \email{yizhen.xu@utah.edu}
#' @source Simulated by Yizhen Xu
"trial_sim_dat"
