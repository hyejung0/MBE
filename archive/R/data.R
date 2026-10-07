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
#' This data set is used in package examples to demonstrate an MBE update
#' without requiring users to refit the historical Stan model.
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



#' Aggregate Leave-One-Out Historical-Model Assessment
#'
#' @description A one-row summary of leave-one-out cross-validation for the
#' random-intercept historical model with inverse-gamma variance priors. Each of
#' the 66 simulated trials was held out once. The tolerance-interval calculation
#' includes both residual clinical heterogeneity and the held-out trial's
#' clinical sampling variance.
#'
#' @details Archived results from the earlier development assessment. RMSE
#' uses full-carryover MBE posterior means, incorporating the held-out clinical
#' estimate. Coverage used empirical-percentile indicators and an unweighted
#' surrogate-update simulation. The current [clinical_predictive_distribution()]
#' uses joint conditioning and surrogate-likelihood weights, and
#' [tolerance_interval_coverage()] checks inclusive quantile bounds. This table
#' is retained unchanged and is not a result of the current coverage function.
#' These are internal cross-validation results on simulated data, not external
#' validation. Model specification here describes only this archived dataset;
#' the current assessment functions do not prescribe an intercept or prior.
#'
#' @format A data frame with 1 row and 4 variables:
#' \describe{
#'   \item{n_trials}{Integer. Number of held-out trials.}
#'   \item{loo_cv_rmse}{Numeric. Root mean squared error between held-out
#'   clinical estimates and their posterior means.}
#'   \item{coverage_95}{Numeric. Proportion of held-out clinical estimates
#'   covered by the 95 percent predictive interval.}
#'   \item{coverage_90}{Numeric. Proportion of held-out clinical estimates
#'   covered by the 90 percent predictive interval.}
#' }
#' @source Archived 66-fit analysis of `trial_sim_dat` with seed 2026, saved in
#' `data-raw/loo-random-inverse-gamma/loo_assessment.rds`. See the preparation
#' script `data-raw/loo_assessment.R` in the source repository.
"loo_assessment_summary"



#' Trial-Level Leave-One-Out Historical-Model Assessment
#'
#' @description Trial-level results underlying `loo_assessment_summary`. The
#' assessment uses the random-intercept historical model with inverse-gamma
#' variance priors and seed 2026.
#'
#' @details These archived results use the earlier empirical-percentile
#' coverage calculation; see [loo_assessment_summary] for provenance and its
#' differences from the current prediction and coverage functions. The paired
#' effects can still be supplied to [loo_cv_rmse()] without refitting any model.
#'
#' @format A data frame with 66 rows and 11 variables:
#' \describe{
#'   \item{trial_id}{Character. Identifier of the held-out trial.}
#'   \item{observed_clinical}{Numeric. Observed clinical treatment-effect estimate.}
#'   \item{posterior_mean_clinical}{Numeric. Posterior mean for the held-out
#'   clinical treatment effect.}
#'   \item{squared_error}{Numeric. Squared error of the posterior mean.}
#'   \item{predictive_percentile}{Numeric. Empirical predictive percentile of
#'   the observed clinical estimate.}
#'   \item{predictive_lower_95}{Numeric. Lower endpoint of the 95 percent
#'   predictive interval.}
#'   \item{predictive_lower_90}{Numeric. Lower endpoint of the 90 percent
#'   predictive interval.}
#'   \item{predictive_upper_90}{Numeric. Upper endpoint of the 90 percent
#'   predictive interval.}
#'   \item{predictive_upper_95}{Numeric. Upper endpoint of the 95 percent
#'   predictive interval.}
#'   \item{covered_95}{Logical. Whether the 95 percent interval covered the
#'   observed clinical estimate.}
#'   \item{covered_90}{Logical. Whether the 90 percent interval covered the
#'   observed clinical estimate.}
#' }
#' @source Archived 66-fit analysis of `trial_sim_dat` with seed 2026, saved in
#' `data-raw/loo-random-inverse-gamma/loo_assessment.rds`. See the preparation
#' script `data-raw/loo_assessment.R` in the source repository.
"loo_assessment_by_trial"
