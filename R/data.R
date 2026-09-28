#' Simulated Data for Interim Analysis
#'
#' @description This is a simulated data set for chronic kidney disease clinical trials.
#' The definitive clinical endpoint is the first of ESRD or a 57% decline from baseline eGFR.
#' In addition to the treatment effect estimated on the definitive clinical endpoint,
#' the treatment effect is estimated on exactly two surrogate endpoints: chronic
#' eGFR slope (`Sur1`) and acute eGFR slope (`Sur2`).
#'
#' The data contain eight scenarios. Cases 1 and 2 use a baseline eGFR range of
#' 30--60, cases 3 and 4 use 30--75, and cases 5--8 use 25--50
#' ml/min/1.73m2. Odd-numbered cases have no slope treatment effect; even-numbered
#' cases have an active slope treatment effect of 0.8. Cases 7 and 8 repeat the
#' settings of cases 5 and 6, respectively.
#'
#' @format A data frame with 528 rows and 11 variables:
#' \describe{
#'   \item{case}{Integer. The scenario case number (1 through 8).}
#'   \item{analysis.month}{Integer. Number of months since the study began for the specific case.}
#'   \item{R1Clin}{Numeric. Estimated correlation between chronic slope and clinical endpoint. If NA, that is because the treatment effect on clinical endpoint (ClnEst and ClnSE) were not estimated due to small event size.}
#'   \item{R2Clin}{Numeric. Estimated correlation between acute slope and clinical endpoint. If NA, that is because the treatment effect on clinical endpoint (ClnEst and ClnSE) were not estimated due to small event size.}
#'   \item{R12}{Numeric. Estimated correlation between chronic slope and acute slope.}
#'   \item{ClnEst}{Numeric. Estimated treatment effect on the clinical endpoint (log hazard ratio). If NA, that is because the treatment effect on clinical endpoint (ClnEst and ClnSE) were not estimated due to small event size.}
#'   \item{ClnSE}{Numeric. Standard error of the estimated treatment effect on the clinical endpoint. If NA, that is because the clinical treatment effect was not estimated due to small event size.}
#'   \item{Sur2Est}{Numeric. Estimated treatment effect on acute slope.}
#'   \item{Sur2SE}{Numeric. Standard error of the estimated treatment effect on acute slope.}
#'   \item{Sur1Est}{Numeric. Estimated treatment effect on chronic slope.}
#'   \item{Sur1SE}{Numeric. Standard error of the estimated treatment effect on chronic slope.}
#' }
#' @author Jian Ying \email{jian.ying@@hsc.utah.edu}
#' @source Simulated by Jian Ying
"interim_sim_dat"





#' Historical Posterior for Two Surrogate Endpoints
#'
#' Posterior draws of all model parameters from the MCMC model with chronic
#' and acute eGFR slopes as the two surrogates for the clinical endpoint.
#' This data set is used in the examples of the package to demonstrate how
#' to use the `vector_to_matrix` and `MBE` functions.
#'
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
