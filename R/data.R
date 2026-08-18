#' Simulated Data for Interim Analysis
#'
#' @description This is a simulated data set for chronic kidney disease clinical trials.
#' The definitive clinical endpoint is the first of ESRD or a 57% decline from baseline eGFR.
#' In addition to the treatment effect estimated on the definitive clinical endpoint,
#' the treatment effect is also estimated on two surrogate endpoints: the acute slope (Sur2) and chronic slope (Sur1) of eGFR, and on log(UACR) (Sur3).
#' It contains 6 different scenarios:
#' Scenarios 2 is where the baseline eGFR ranges from 30-60 ml/min/1.73m2 (standard baseline eGFR range),
#' scenarios 4 is where the baseline eGFR ranges from 30-70 ml/min/1.73m2 (wider/higher baseline eGFR range),
#' and scenarios 6 is where the baseline eGFR ranges from 25-50 ml/min/1.73m2 (lower baseline eGFR range).
#' For all 3 scenarios, the slope effect is 0.8 and treatment mean log(UACR) change is -0.3.
#' Thus, in all three cases, there is an active treatment effect.
#' Scenario 1 is same as scenario 2, but with no treatment effect (slope effect = 0 and mean log(UACR) change = 0).
#' Likewise, scenario 3 is same as scenario 4, but with no treatment effect, and scenario 5 is same as scenario 6, but with no treatment effect.
#'
#' @format A data frame with 63 rows and 16 variables:
#' \describe{
#'   \item{case}{Integer. The scenario case number (2, 4, or 6).}
#'   \item{analysis.month}{Integer. Number of months since the study began for the specific case.}
#'   \item{R1Clin}{Numeric. Estimated correlation between chronic slope and clinical endpoint. If NA, that is because the treatment effect on clinical endpoint (ClnEst and ClnSE) were not estimated due to small event size.}
#'   \item{R2Clin}{Numeric. Estimated correlation between acute slope and clinical endpoint. If NA, that is because the treatment effect on clinical endpoint (ClnEst and ClnSE) were not estimated due to small event size.}
#'   \item{R3Clin}{Numeric. Estimated correlation between log UACR and clinical endpoint. If NA, that is because the treatment effect on clinical endpoint (ClnEst and ClnSE) were not estimated due to small event size.}
#'   \item{R12}{Numeric. Estimated correlation between chronic slope and acute slope.}
#'   \item{R13}{Numeric. Estimated correlation between chronic slope and log UACR.}
#'   \item{R23}{Numeric. Estimated correlation between acute slope and log UACR.}
#'   \item{ClnEst}{Numeric. Estimated treatment effect on the clinical endpoint (log hazard ratio). If NA, that is because the treatment effect on clinical endpoint (ClnEst and ClnSE) were not estimated due to small event size.}
#'   \item{ClnSE}{Numeric. Standard error of the estimated treatment effect on the clinical endpoint.. If NA, that is because the treatment effect on clinical endpoint (ClnEst and ClnSE) were not estimated due to small event size.}
#'   \item{Sur3Est}{Numeric. Estimated treatment effect on log UACR.}
#'   \item{Sur3SE}{Numeric. Standard error of the estimated treatment effect on log UACR.}
#'   \item{Sur2Est}{Numeric. Estimated treatment effect on acute slope.}
#'   \item{Sur2SE}{Numeric. Standard error of the estimated treatment effect on acute slope.}
#'   \item{Sur1Est}{Numeric. Estimated treatment effect on chronic slope.}
#'   \item{Sur1SE}{Numeric. Standard error of the estimated treatment effect on chronic slope.}
#' }
#' @author Jian Ying \email{jian.ying@@hsc.utah.edu}
#' @source Simulated by Jian Ying
"interim_sim_dat"





#' Posterior Distribution of Parameter Values from Modeling with UACR, Acute Slope, and Chronic Slope as Surrogates for Clinical Endpoint
#'
#' Posterior distribution of all parameters from the MCMC model with UACR,
#' acute slope, and chronic slope as surrogates for the clinical endpoint.
#' This data set is used in the examples of the package to demonstrate how
#' to use the `vector_to_matrix` and `MBE` functions.
#'
#' @format A data.table with 500 rows and 4 variables:
#' \describe{
#'   \item{col1}{Description}
#'   \item{col2}{Description}
#' }
"historical_posterior"




#' 66 Trial Level Summary Simulated Data
#'
#' @description Wait until Yizhen gives me her description
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
