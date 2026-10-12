#' Predict an Observed Clinical Effect from Two Observed Surrogate Effects
#'
#' @description Generates equally weighted predictive draws for a new
#' trial's observed clinical treatment-effect estimate using historical
#' posterior and observed treatment effects on the surrogate endpoints.
#'
#' @param mcmc_dat Historical posterior draws as `posterior_samples` returned by
#' [historical_model_fit_2surrogates()] function.
#' @param sample_dat A named list containing `ClnSE`, `Sur1Est`, `Sur1SE`,
#'   `Sur2Est`, `Sur2SE`, `R1Clin`, `R2Clin`, and `R12` for the one trial for which
#'   we wish to estimate the posterior predictive distribution.
#' @param seed Optional positive integer for reproducible simulation. When
#'   supplied, the caller's random-number state is restored on exit.
#' @param diffuse Logical; replace the historical surrogate-effect distribution
#'   with a diffuse prior, as in [surrogate_MBE()]. Defaults to `FALSE` (full carryover).
#' @param diffuse_se Positive standard deviation of the diffuse surrogate prior.
#' @param intercept0 Logical; center the clinical-intercept draws at zero, as in
#'   [surrogate_MBE()]. Defaults to `FALSE`. An absent clinical-intercept column in a
#'   fixed-intercept fit is treated as zero.
#'
#' @details Review the details provided in [surrogate_MBE()]. We can use the
#' historical posterior draws and the observed treatment effect on the surrogate
#' effects to sample posterior distribution on the true treatment effect on the
#' clinical endpoint. Note also
#'\deqn{\hat{\theta}_0 \sim N(\theta_0, \sigma^2_0)} where \eqn{\sigma^2_0} is
#' the known sampling variance of the clinical treatment effect
#' (i.e., `ClnSE` of the new trial). We can then collect the predictive
#' distribution on the observed clinical treatment effect by sampling the error
#' term from \eqn{N(0, \sigma^2_0)} and adding it to the surrogate MBE's posterior
#' distribution of the true clinical effect.
#'
#'
#' @return A numeric vector of posterior predictive distribution of the treatment
#' effect on the clinical endpoint. The length is same as the number of rows of
#' the historical posteriors.
#' @export
#' @examples
#' data("historical_posterior")
#' data("trial_sim_dat")
#' held_out <- list(
#'   ClnSE = trial_sim_dat$CE_se[1],
#'   Sur1Est = trial_sim_dat$Sur1_est[1],
#'   Sur1SE = trial_sim_dat$Sur1_se[1],
#'   Sur2Est = trial_sim_dat$Sur2_est[1],
#'   Sur2SE = trial_sim_dat$Sur2_se[1],
#'   R1Clin = trial_sim_dat$Cor_CE_Sur1[1],
#'   R2Clin = trial_sim_dat$Cor_CE_Sur2[1],
#'   R12 = trial_sim_dat$Cor_Sur1_Sur2[1]
#' )
#' predictive <- clinical_predictive_distribution(
#'   mcmc_dat=historical_posterior,
#'   sample_dat=held_out,
#'   seed = 2026
#' )
#' stats::quantile(predictive, c(0.025, 0.975))

clinical_predictive_distribution <- function(mcmc_dat, sample_dat,
                                            seed = NULL,
                                             diffuse = FALSE, diffuse_se = 100,
                                             intercept0 = FALSE) {

  #First, get the posterior distribution of the true treatment effect on the
  #clinical endpoint using surrogate_MBE

  MBE_from_surrogate_only<-surrogate_MBE(
   mcmc_dat = mcmc_dat,
   sample_dat=sample_dat,
   diffuse = diffuse,
   diffuse_se = diffuse_se,
   intercept0 = intercept0)

  #generate random error using the standard error of the observed
  #treatment effect on the clinical endpoint
  error <- stats::rnorm(nrow(mcmc_dat), 0, sample_dat$ClnSE)

  # Add the random error and return it
  return(
    MBE_from_surrogate_only$post_CE + error
  )
}
