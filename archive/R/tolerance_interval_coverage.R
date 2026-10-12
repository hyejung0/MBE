#' Check Capture of an Observed Clinical Effect in Predictive Intervals
#'
#' @description Calls [clinical_predictive_distribution()] using an already
#' fitted historical posterior, then checks whether the observed clinical
#' estimate falls within user-specified predictive intervals. Historical models
#' are not refitted. The observed clinical estimate is used only for capture.
#'
#' @param observed One finite numeric observed clinical treatment-effect estimate.
#' @inheritParams clinical_predictive_distribution
#' @param probs Lower and upper quantile probabilities: either a numeric vector
#'   of length two, such as `c(0.025, 0.975)`, or a numeric two-column matrix
#'   (or data frame) with one row per interval. Column 1 is the lower probability
#'   and column 2 the upper probability. Each pair must satisfy
#'   `0 <= lower < upper <= 1`. Asymmetric intervals are allowed.
#' @param ... Additional arguments to [clinical_predictive_distribution()],
#'   including `n_draws`, `seed`, `diffuse`, `diffuse_se`, and `intercept0`.
#'
#' @details Predictive draws are generated once and shared across all intervals
#' in a call. Bounds use [stats::quantile()] with `type = 7`. Capture is inclusive:
#' `lower <= observed && observed <= upper`. Missing or infinite inputs are
#' rejected. Results follow the supplied row order. Specify `seed` for
#' reproducibility; finite Monte Carlo draws can affect near-boundary results.
#'
#' The name follows the manuscript's tolerance-interval terminology. These are
#' posterior predictive intervals, not classical tolerance intervals with a
#' separate confidence level. Repeat this function across held-out trials and
#' average `covered` for each probability pair to obtain aggregate coverage.
#'
#' Earlier development code used empirical percentiles and a different
#' predictive simulation. Archived `loo_assessment_by_trial` indicators retain
#' that earlier calculation; see [clinical_predictive_distribution()] for the
#' current conditioning and weighting method.
#'
#' @return A data frame with one row per interval: `lower_prob`, `upper_prob`,
#'   `level` (100 times their difference), `observed`, interval bounds `lower`
#'   and `upper`, and logical `covered`. For one trial, `covered` is a capture
#'   indicator rather than an aggregate coverage rate.
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
#' tolerance_interval_coverage(
#'   observed = trial_sim_dat$CE_est[1],
#'   mcmc_dat = historical_posterior[1:100, ], sample_dat = held_out,
#'   probs = rbind(c(0.025, 0.975), c(0.05, 0.95), c(0.075, 0.925)),
#'   seed = 2026
#' )
tolerance_interval_coverage <- function(observed, mcmc_dat, sample_dat,
                                        probs = c(0.025, 0.975), ...) {
  if (!is.numeric(observed) || length(observed) != 1L ||
      !is.null(dim(observed)) || !is.finite(observed)) {
    stop("`observed` must be one finite numeric value.", call. = FALSE)
  }
  bounds <- .validate_interval_probs(probs)
  draws <- clinical_predictive_distribution(mcmc_dat, sample_dat, ...)
  .predictive_interval_capture(observed, draws, bounds)
}

.validate_interval_probs <- function(probs) {
  if (is.numeric(probs) && is.null(dim(probs)) && length(probs) == 2L) {
    probs <- matrix(probs, nrow = 1L)
  }
  if (is.data.frame(probs)) probs <- as.matrix(probs)
  if (!is.matrix(probs) || !is.numeric(probs) || ncol(probs) != 2L ||
      nrow(probs) < 1L || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1) || any(probs[, 1] >= probs[, 2])) {
    stop("`probs` must supply lower/upper probability pairs with 0 <= lower < upper <= 1.",
         call. = FALSE)
  }
  probs
}

.predictive_interval_capture <- function(observed, draws, probs) {
  lower <- stats::quantile(draws, probs[, 1], names = FALSE, type = 7)
  upper <- stats::quantile(draws, probs[, 2], names = FALSE, type = 7)
  data.frame(
    lower_prob = unname(probs[, 1]), upper_prob = unname(probs[, 2]),
    level = unname(100 * (probs[, 2] - probs[, 1])),
    observed = as.numeric(observed), lower = lower, upper = upper,
    covered = observed >= lower & observed <= upper
  )
}
