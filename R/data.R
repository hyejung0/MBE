#' Simulated Data for Package Examples
#'
#' @description This is a simulated data set for chronic kidney disease clinical trials.
#' The definitive clinical endpoint is the first of ESRD or a 57% decline from baseline eGFR.
#' It has treatment effect on chronic slope(Sur1), acute slope (Sur2), and UACR (Sur3),
#' for intermediate analysis months 15-35. It contains 3 different scenarios
#' Scenarios 2 is where the baseline eGFR ranges from 30-60 ml/min/1.73m2 (standard baselineeGFR range),
#' scenarios 4 is where the baseline eGFR ranges from 30-70 ml/min/1.73m2 (wider/higher baseline eGFR range),
#' and scenarios 6 is where the baseline eGFR ranges from 25-50 ml/min/1.73m2 (lower baseline eGFR range).
#' For all 3 scenarios, the slope effect is 0.8 and treatment mean log(UACR) change is -0.3.
#' Thus, in all three cases, there is an active treatment effect.
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
"sim_dat"

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
