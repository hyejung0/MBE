#' Predict an Observed Clinical Effect from Two Observed Surrogate Effects
#'
#' @description Generates equally weighted predictive draws for a new
#' trial's observed clinical treatment-effect estimate using an already fitted
#' historical posterior. This function does not fit a historical model and does
#' not use the held-out clinical estimate to construct predictions.
#'
#' @param mcmc_dat Historical posterior draws in the format accepted by [MBE()].
#'   For leave-one-trial-out assessment, fit these draws without the trial to be
#'   predicted, using [historical_model_fit_2surrogates()]. Fixed- and
#'   random-intercept fits and both uncertainty priors are supported.
#' @param sample_dat A named list containing `ClnSE`, `Sur1Est`, `Sur1SE`,
#'   `Sur2Est`, `Sur2SE`, `R1Clin`, `R2Clin`, and `R12` for one held-out trial.
#'   An optional `ClnEst` entry is ignored. These are the same sampling inputs
#'   used by [MBE()], except the clinical estimate is not needed.
#' @param n_draws Number of predictive draws to generate. Defaults to the number
#'   of rows in `mcmc_dat`.
#' @param seed Optional positive integer for reproducible simulation. When
#'   supplied, the caller's random-number state is restored on exit.
#' @param diffuse Logical; replace the historical surrogate-effect distribution
#'   with a diffuse prior, as in [MBE()]. Defaults to `FALSE` (full carryover).
#' @param diffuse_se Positive standard deviation of the diffuse surrogate prior.
#' @param intercept0 Logical; center the clinical-intercept draws at zero, as in
#'   [MBE()]. Defaults to `FALSE`. An absent clinical-intercept column in a
#'   fixed-intercept fit is treated as zero.
#'
#' @details For each historical draw, the joint distribution of the observed
#' clinical and surrogate estimates is multivariate normal, with the model
#' covariance plus the within-trial sampling covariance. The function conditions
#' this joint distribution on the two observed surrogate estimates using
#' normal conditional-distribution formulas to generate predictive draws for the new trial.
#' Historical parameter draws are reweighted by the marginal likelihood of
#' those surrogate estimates, then sampled to form the predictive mixture.
#'
#' This includes residual clinical heterogeneity, clinical sampling variance,
#' surrogate measurement uncertainty, posterior parameter uncertainty, and the
#' supplied within-trial sampling correlations. The clinical estimate itself
#' from the new trial is never used for weighting or conditioning.
#' The returned draws target an observed clinical estimate, not a latent true
#' clinical effect.
#'
#' This joint conditioning also accounts for clinical-surrogate sampling
#' correlations and updates historical-draw weights. It differs from the earlier
#' development assessment's unweighted surrogate-update simulation, so new
#' predictive intervals need not reproduce archived coverage results.
#'
#' @return A numeric vector of `n_draws` equally weighted predictive draws.
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
#'   historical_posterior[1:100, ], held_out, seed = 2026
#' )
#' stats::quantile(predictive, c(0.025, 0.975))
clinical_predictive_distribution <- function(mcmc_dat, sample_dat,
                                             n_draws = NULL, seed = NULL,
                                             diffuse = FALSE, diffuse_se = 100,
                                             intercept0 = FALSE) {
  .validate_mcmc_dat(mcmc_dat)
  if (!is.list(sample_dat) || is.null(names(sample_dat)) ||
      anyDuplicated(names(sample_dat))) {
    stop("`sample_dat` must be a named list with unique sampling-input names.",
         call. = FALSE)
  }
  # Ignore the clinical observation, even if supplied with an MBE input list.
  prediction_data <- sample_dat
  prediction_data$ClnEst <- 0
  .validate_sample_dat(prediction_data)
  if (is.null(n_draws)) n_draws <- nrow(mcmc_dat)
  .validate_prediction_integer(n_draws, "n_draws")
  if (!is.null(seed)) .validate_prediction_integer(seed, "seed")

  prior <- vector_to_matrix(
    mcmc_dat, diffuse = diffuse, diffuse_se = diffuse_se,
    intercept0 = intercept0
  )
  se <- c(prediction_data$ClnSE, prediction_data$Sur1SE, prediction_data$Sur2SE)
  correlation <- matrix(c(
    1, prediction_data$R1Clin, prediction_data$R2Clin,
    prediction_data$R1Clin, 1, prediction_data$R12,
    prediction_data$R2Clin, prediction_data$R12, 1
  ), nrow = 3)
  sampling_covariance <- correlation * outer(se, se)
  .validate_positive_definite(sampling_covariance, "The sampling covariance matrix")
  observed_surrogates <- c(prediction_data$Sur1Est, prediction_data$Sur2Est)
  components <- .clinical_predictive_components(prior, sampling_covariance, observed_surrogates)

  simulate <- function() {
    index <- sample.int(nrow(mcmc_dat), n_draws, replace = TRUE,
                        prob = components$weight)
    stats::rnorm(n_draws, components$mean[index], sqrt(components$variance[index]))
  }
  .with_preserved_seed(seed, simulate())
}

.clinical_predictive_components <- function(prior, sampling_covariance, observed_surrogates) {
  n <- length(prior$mean)
  mean <- variance <- log_weight <- numeric(n)
  for (i in seq_len(n)) {
    mu <- as.numeric(prior$mean[[i]])
    covariance <- prior$covar[[i]] + sampling_covariance
    surrogate_covariance <- covariance[2:3, 2:3, drop = FALSE]
    regression <- covariance[1, 2:3, drop = FALSE] %*%
      chol2inv(chol(surrogate_covariance))
    mean[i] <- mu[1] + as.numeric(regression %*% (observed_surrogates - mu[2:3]))
    variance[i] <- covariance[1, 1] - as.numeric(regression %*% covariance[2:3, 1])
    log_weight[i] <- mvtnorm::dmvnorm(
      observed_surrogates, mu[2:3], surrogate_covariance, log = TRUE
    )
  }
  if (any(!is.finite(mean)) || any(!is.finite(variance)) || any(variance <= 0) ||
      anyNA(log_weight) || any(log_weight == Inf) || !any(is.finite(log_weight))) {
    stop("Cannot construct a finite predictive distribution from these inputs.",
         call. = FALSE)
  }
  weight <- exp(log_weight - max(log_weight))
  list(mean = mean, variance = variance, weight = weight / sum(weight))
}

.validate_prediction_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x > .Machine$integer.max || x != floor(x)) {
    stop(sprintf("`%s` must be one positive integer within R's integer range.", name),
         call. = FALSE)
  }
}

.with_preserved_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}
