.mbe_loo_draw_fields <- c(
  "alphaCEonSur1Sur2",
  "b1CEonSur1Sur2",
  "b2CEonSur1Sur2",
  "SigSqCEonSur1Sur2",
  "alphaSur1onSur2",
  "bSur1onSur2",
  "SigSqSur1onSur2",
  "muSur2",
  "sigSqSur2"
)

.prepare_loo_data <- function(data, trial_id = NULL, minimum_rows = 1L) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.", call. = FALSE)
  }

  data <- as.data.frame(data)
  if (nrow(data) < minimum_rows) {
    stop(
      sprintf("`data` must contain at least %d trials.", minimum_rows),
      call. = FALSE
    )
  }

  missing_fields <- setdiff(.mbe_historical_fields, names(data))
  if (length(missing_fields)) {
    stop(
      paste0(
        "`data` is missing fields required for the two-surrogate model: ",
        .format_field_list(missing_fields), "."
      ),
      call. = FALSE
    )
  }

  historical_data <- data[, .mbe_historical_fields, drop = FALSE]
  .validate_historical_data(historical_data)

  if (is.null(trial_id)) {
    trial_id <- if ("trial_id" %in% names(data)) "trial_id" else NULL
  }
  if (!is.null(trial_id)) {
    if (!is.character(trial_id) || length(trial_id) != 1L ||
        is.na(trial_id) || !nzchar(trial_id)) {
      stop("`trial_id` must be NULL or one column name.", call. = FALSE)
    }
    if (!trial_id %in% names(data)) {
      stop(sprintf("Trial ID column `%s` was not found in `data`.", trial_id), call. = FALSE)
    }
    ids <- data[[trial_id]]
  } else {
    ids <- seq_len(nrow(data))
  }

  if (anyNA(ids) || anyDuplicated(ids)) {
    stop("Trial identifiers must be non-missing and unique.", call. = FALSE)
  }

  list(
    data = data,
    historical_data = historical_data,
    ids = as.character(ids)
  )
}

.historical_row_to_sample_dat <- function(data, i) {
  list(
    ClnEst = data$CE_est[i],
    ClnSE = data$CE_se[i],
    Sur1Est = data$Sur1_est[i],
    Sur1SE = data$Sur1_se[i],
    Sur2Est = data$Sur2_est[i],
    Sur2SE = data$Sur2_se[i],
    R1Clin = data$Cor_CE_Sur1[i],
    R2Clin = data$Cor_CE_Sur2[i],
    R12 = data$Cor_Sur1_Sur2[i]
  )
}

.with_preserved_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  .validate_positive_integer(seed, "seed")

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
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

.clinical_predictive_sd <- function(residual_variance, clinical_se) {
  sqrt(residual_variance + clinical_se^2)
}

.extract_loo_draws <- function(x, expected_id = NULL) {
  source_label <- "the supplied object"
  if (is.character(x) && length(x) == 1L && !is.na(x)) {
    source_label <- sprintf("`%s`", x)
    if (!file.exists(x)) {
      stop(sprintf("LOO posterior file %s does not exist.", source_label), call. = FALSE)
    }
    x <- readRDS(x)
  }

  if (is.list(x) && !is.data.frame(x) && !is.matrix(x) &&
      !is.null(x$posterior_draws)) {
    if (!is.null(expected_id) && !is.null(x$held_out_trial) &&
        !identical(as.character(x$held_out_trial), as.character(expected_id))) {
      stop(
        sprintf(
          "LOO posterior for trial `%s` was supplied for trial `%s`.",
          x$held_out_trial, expected_id
        ),
        call. = FALSE
      )
    }
    x <- x$posterior_draws
  }

  if (!is.data.frame(x) && !is.matrix(x)) {
    stop(
      sprintf("LOO posterior %s must contain a data frame or matrix of draws.", source_label),
      call. = FALSE
    )
  }
  .validate_mcmc_dat(x)
  x
}

.validate_saved_loo_record <- function(path, expected_id, expected_seed) {
  record <- tryCatch(
    readRDS(path),
    error = function(e) {
      stop(
        sprintf(
          "Existing LOO posterior `%s` could not be read; remove it or use `overwrite = TRUE`.",
          path
        ),
        call. = FALSE
      )
    }
  )
  expected_specification <- list(
    random_intercept = TRUE,
    prior_for_uncertainty = "inverse_gamma"
  )

  if (!is.list(record) ||
      !identical(as.character(record$held_out_trial), as.character(expected_id)) ||
      !isTRUE(all.equal(as.numeric(record$seed), as.numeric(expected_seed))) ||
      !identical(record$specification, expected_specification) ||
      is.null(record$posterior_draws) ||
      length(setdiff(.mbe_loo_draw_fields, names(record$posterior_draws)))) {
    stop(
      sprintf(
        paste0(
          "Existing LOO posterior `%s` does not match the requested trial, seed, ",
          "or model specification; remove it or use `overwrite = TRUE`."
        ),
        path
      ),
      call. = FALSE
    )
  }
  .extract_loo_draws(record, expected_id = expected_id)
  invisible(record)
}

.loo_trial_metrics <- function(mcmc_dat, sample_dat) {
  draws <- as.data.frame(mcmc_dat)
  if (!"alphaCEonSur1Sur2" %in% names(draws)) {
    draws$alphaCEonSur1Sur2 <- 0
  }

  # Match the manuscript's full-carryover MBE calculation for LOO-CV RMSE.
  mbe_result <- MBE(
    mcmc_dat = draws,
    sample_dat = sample_dat,
    diffuse = FALSE,
    intercept0 = FALSE
  )
  posterior_mean_clinical <- as.numeric(mbe_result$post_mean[1, 1])
  squared_error <- (sample_dat$ClnEst - posterior_mean_clinical)^2

  # For tolerance coverage, condition only on the two observed surrogates.
  surrogate_observed <- matrix(c(sample_dat$Sur1Est, sample_dat$Sur2Est), ncol = 1)
  surrogate_sampling_covar <- matrix(
    c(
      sample_dat$Sur1SE^2,
      sample_dat$R12 * sample_dat$Sur1SE * sample_dat$Sur2SE,
      sample_dat$R12 * sample_dat$Sur1SE * sample_dat$Sur2SE,
      sample_dat$Sur2SE^2
    ),
    nrow = 2,
    byrow = TRUE
  )
  .validate_positive_definite(
    surrogate_sampling_covar,
    "The surrogate sampling covariance matrix"
  )
  surrogate_sampling_precision <- chol2inv(chol(surrogate_sampling_covar))

  prior <- vector_to_matrix(
    MCMC_dat = draws,
    diffuse = FALSE,
    intercept0 = FALSE
  )
  B <- nrow(draws)
  predictive_clinical <- numeric(B)

  for (b in seq_len(B)) {
    surrogate_prior_mean <- prior$mean[[b]][2:3, , drop = FALSE]
    surrogate_prior_covar <- prior$covar[[b]][2:3, 2:3, drop = FALSE]
    surrogate_prior_precision <- chol2inv(chol(surrogate_prior_covar))

    surrogate_post_covar <- chol2inv(chol(
      surrogate_prior_precision + surrogate_sampling_precision
    ))
    surrogate_post_mean <- surrogate_post_covar %*% (
      surrogate_prior_precision %*% surrogate_prior_mean +
        surrogate_sampling_precision %*% surrogate_observed
    )
    surrogate_draw <- as.numeric(mvtnorm::rmvnorm(
      1,
      mean = surrogate_post_mean,
      sigma = surrogate_post_covar
    ))

    clinical_mean <-
      draws$alphaCEonSur1Sur2[b] +
      draws$b1CEonSur1Sur2[b] * surrogate_draw[1] +
      draws$b2CEonSur1Sur2[b] * surrogate_draw[2]

    # Corrected from the original analysis script: prediction of the observed
    # clinical estimate includes both between-trial residual heterogeneity and
    # the held-out trial's clinical sampling variance.
    predictive_clinical[b] <- stats::rnorm(
      1,
      mean = clinical_mean,
      sd = .clinical_predictive_sd(
        draws$SigSqCEonSur1Sur2[b],
        sample_dat$ClnSE
      )
    )
  }

  predictive_percentile <- mean(predictive_clinical <= sample_dat$ClnEst)
  predictive_quantiles <- stats::quantile(
    predictive_clinical,
    probs = c(0.025, 0.05, 0.95, 0.975),
    names = FALSE
  )

  list(
    observed_clinical = sample_dat$ClnEst,
    posterior_mean_clinical = posterior_mean_clinical,
    squared_error = squared_error,
    predictive_percentile = predictive_percentile,
    predictive_lower_95 = predictive_quantiles[1],
    predictive_lower_90 = predictive_quantiles[2],
    predictive_upper_90 = predictive_quantiles[3],
    predictive_upper_95 = predictive_quantiles[4],
    covered_95 = predictive_percentile >= 0.025 && predictive_percentile <= 0.975,
    covered_90 = predictive_percentile >= 0.05 && predictive_percentile <= 0.95
  )
}

#' Fit Random-Intercept Inverse-Gamma Models for Leave-One-Out Assessment
#'
#' @description Fits the package's initial manuscript specification once for
#' each held-out trial: a random clinical intercept with inverse-gamma priors
#' for variance parameters. Each fit uses all remaining trials. Only the nine
#' posterior columns needed for MBE assessment are saved, rather than complete
#' CmdStan fit objects or chain CSV files.
#'
#' @param data A data frame or matrix containing the nine historical-model
#'   fields documented in [historical_model_fit_2surrogates()]. Additional
#'   columns, such as `trial_id` and simulation truths, are allowed.
#' @param output_dir Directory in which to save one compact `.rds` file per
#'   held-out trial. A directory under `data-raw/` is recommended; the full
#'   files should not be included in the installed package.
#' @param trial_id Optional name of a column containing unique trial IDs. When
#'   `NULL`, `trial_id` is used if present; otherwise row numbers are used.
#' @param seed Positive integer base seed. Trial `i` uses `seed + i - 1`.
#' @param overwrite Logical; overwrite existing trial files when `TRUE`.
#' @param nchains Number of MCMC chains for each historical fit.
#' @param ncores Number of chains to run in parallel within each fit.
#' @param niter Number of post-warmup draws per chain.
#' @param nwarmup Number of warmup iterations per chain.
#' @param show_messages Logical; show CmdStan sampling messages.
#' @param ... Additional arguments forwarded to
#'   [historical_model_fit_2surrogates()] and then to `CmdStanModel$sample()`,
#'   such as `adapt_delta`, `max_treedepth`, and `refresh`.
#'
#' @return Invisibly, a named character vector of saved `.rds` paths. Names
#'   are the held-out trial identifiers. Existing files are returned without
#'   refitting when `overwrite = FALSE`, after their trial identifier, seed,
#'   model specification, and required posterior columns are validated.
#' @export
#'
#' @examples
#' \dontrun{
#' paths <- fit_loo_historical_models(
#'   trial_sim_dat,
#'   output_dir = "data-raw/loo-random-inverse-gamma",
#'   seed = 2026,
#'   nchains = 4,
#'   ncores = 4,
#'   niter = 2000,
#'   nwarmup = 1000
#' )
#' }
fit_loo_historical_models <- function(data,
                                      output_dir,
                                      trial_id = NULL,
                                      seed = 1L,
                                      overwrite = FALSE,
                                      nchains = 4L,
                                      ncores = 1L,
                                      niter = 2000L,
                                      nwarmup = 1000L,
                                      show_messages = TRUE,
                                      ...) {
  prepared <- .prepare_loo_data(data, trial_id, minimum_rows = 2L)
  .validate_positive_integer(seed, "seed")
  .validate_scalar_logical(overwrite, "overwrite")
  .validate_scalar_logical(show_messages, "show_messages")
  .validate_positive_integer(nchains, "nchains")
  .validate_positive_integer(ncores, "ncores")
  .validate_positive_integer(niter, "niter")
  .validate_positive_integer(nwarmup, "nwarmup")

  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("`output_dir` must be one non-empty directory path.", call. = FALSE)
  }
  if (file.exists(output_dir) && !dir.exists(output_dir)) {
    stop("`output_dir` exists but is not a directory.", call. = FALSE)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  safe_ids <- gsub("[^A-Za-z0-9._-]+", "_", prepared$ids)
  if (anyDuplicated(safe_ids)) {
    stop("Trial IDs do not produce unique safe filenames.", call. = FALSE)
  }
  paths <- file.path(output_dir, paste0("heldout_trial_", safe_ids, ".rds"))
  names(paths) <- prepared$ids

  for (i in seq_len(nrow(prepared$historical_data))) {
    trial_seed <- seed + i - 1L
    if (trial_seed > .Machine$integer.max) {
      stop("`seed + number of trials` exceeds the supported integer range.", call. = FALSE)
    }
    if (file.exists(paths[i]) && !overwrite) {
      .validate_saved_loo_record(paths[i], prepared$ids[i], trial_seed)
      if (show_messages) {
        message(sprintf("Using existing LOO posterior for trial `%s`.", prepared$ids[i]))
      }
      next
    }

    if (show_messages) {
      message(sprintf(
        "Fitting random-intercept inverse-gamma model with trial `%s` held out (%d/%d).",
        prepared$ids[i], i, nrow(prepared$historical_data)
      ))
    }

    fit <- historical_model_fit_2surrogates(
      data = prepared$historical_data[-i, , drop = FALSE],
      random_intercept = TRUE,
      prior_for_uncertainty = "inverse_gamma",
      output_dir = tempdir(),
      show_messages = show_messages,
      nchains = nchains,
      ncores = ncores,
      niter = niter,
      nwarmup = nwarmup,
      seed = trial_seed,
      ...
    )

    posterior_draws <- as.data.frame(fit$posterior_samples)
    missing_draws <- setdiff(.mbe_loo_draw_fields, names(posterior_draws))
    if (length(missing_draws)) {
      stop(
        paste0(
          "The fitted model did not return required posterior fields: ",
          .format_field_list(missing_draws), "."
        ),
        call. = FALSE
      )
    }
    posterior_draws <- posterior_draws[, .mbe_loo_draw_fields, drop = FALSE]

    saveRDS(
      list(
        held_out_trial = prepared$ids[i],
        specification = list(
          random_intercept = TRUE,
          prior_for_uncertainty = "inverse_gamma"
        ),
        seed = trial_seed,
        posterior_draws = posterior_draws
      ),
      file = paths[i],
      compress = "xz"
    )
  }

  invisible(paths)
}

#' Calculate LOO-CV RMSE and Tolerance-Interval Coverage
#'
#' @description Calculates the historical-model assessment measures described
#' in the MBE manuscript from posterior draws obtained after holding out each
#' trial. LOO-CV RMSE uses the posterior clinical mean from a full-carryover
#' MBE update. Tolerance coverage predicts the held-out observed clinical
#' estimate using its two observed surrogate estimates without conditioning on
#' the observed clinical estimate itself.
#'
#' @param loo_posteriors A list with one posterior-draw data frame or matrix per
#'   held-out trial, or a character vector of `.rds` paths created by
#'   [fit_loo_historical_models()]. The output directory itself may also be
#'   supplied, in which case its saved trial files are discovered and matched
#'   using their metadata. A single data frame or matrix is accepted when
#'   `data` contains one trial. Named inputs are matched to trial IDs; unnamed
#'   inputs must follow the row order of `data`.
#' @param data A data frame or matrix containing the held-out trial observations.
#'   It must include the nine historical-model fields; additional columns are
#'   allowed.
#' @param trial_id Optional name of the unique trial-ID column. When `NULL`,
#'   `trial_id` is used if present; otherwise row numbers are used.
#' @param seed Optional positive integer used for reproducible posterior
#'   predictive simulation. The caller's random-number state is preserved.
#'
#' @details For each MCMC draw, the two true surrogate effects are updated using
#' the two held-out surrogate estimates. A predictive draw for the observed
#' clinical estimate is then generated with variance
#' `SigSqCEonSur1Sur2 + ClnSE^2`. This corrects the original analysis script,
#' which inadvertently omitted `SigSqCEonSur1Sur2` from the predictive variance.
#'
#' The 90% and 95% coverage indicators follow the manuscript calculation: the
#' empirical percentile of the observed clinical estimate must fall within
#' 0.05--0.95 or 0.025--0.975, respectively.
#'
#' @return A list with two data frames:
#' \describe{
#'   \item{per_trial}{Held-out observations, posterior clinical means, squared
#'   errors, predictive percentiles and intervals, and coverage indicators.}
#'   \item{summary}{Number of trials, LOO-CV RMSE, and aggregate 90% and 95%
#'   tolerance-interval coverage.}
#' }
#' @export
#'
#' @examples
#' data("historical_posterior")
#' data("trial_sim_dat")
#'
#' # The packaged historical posterior was fit with trial 1 held out, so it can
#' # demonstrate the calculation for one trial. Aggregate coverage requires one
#' # separately fitted posterior for every held-out trial.
#' one_trial <- loo_cv_model_assessment(
#'   loo_posteriors = list(`1` = historical_posterior[1:100, ]),
#'   data = trial_sim_dat[1, ],
#'   seed = 1
#' )
#' one_trial$per_trial
#' one_trial$summary
loo_cv_model_assessment <- function(loo_posteriors,
                                    data,
                                    trial_id = NULL,
                                    seed = NULL) {
  prepared <- .prepare_loo_data(data, trial_id, minimum_rows = 1L)

  if (is.character(loo_posteriors) && length(loo_posteriors) == 1L &&
      !is.na(loo_posteriors) && dir.exists(loo_posteriors)) {
    posterior_files <- list.files(
      loo_posteriors,
      pattern = "^heldout_trial_.*[.]rds$",
      full.names = TRUE
    )
    if (!length(posterior_files)) {
      stop("The LOO posterior directory contains no held-out trial files.", call. = FALSE)
    }
    posterior_ids <- vapply(posterior_files, function(path) {
      record <- readRDS(path)
      if (!is.list(record) || is.null(record$held_out_trial)) {
        stop(
          sprintf("LOO posterior file `%s` has no held-out trial metadata.", path),
          call. = FALSE
        )
      }
      as.character(record$held_out_trial)
    }, character(1))
    names(posterior_files) <- posterior_ids
    loo_posteriors <- posterior_files
  }

  if (is.data.frame(loo_posteriors) || is.matrix(loo_posteriors)) {
    loo_posteriors <- list(loo_posteriors)
  } else if (is.character(loo_posteriors)) {
    input_names <- names(loo_posteriors)
    loo_posteriors <- as.list(loo_posteriors)
    names(loo_posteriors) <- input_names
  }
  if (!is.list(loo_posteriors)) {
    stop(
      "`loo_posteriors` must be posterior draws, a list of draws, or saved `.rds` paths.",
      call. = FALSE
    )
  }
  if (length(loo_posteriors) != nrow(prepared$historical_data)) {
    stop(
      sprintf(
        "`loo_posteriors` has %d elements but `data` has %d trials.",
        length(loo_posteriors), nrow(prepared$historical_data)
      ),
      call. = FALSE
    )
  }

  posterior_names <- names(loo_posteriors)
  if (!is.null(posterior_names)) {
    if (any(!nzchar(posterior_names)) || anyDuplicated(posterior_names) ||
        !setequal(posterior_names, prepared$ids)) {
      stop(
        "Named `loo_posteriors` must have one unique name matching every trial ID.",
        call. = FALSE
      )
    }
    loo_posteriors <- loo_posteriors[prepared$ids]
  }

  calculate <- function() {
    rows <- lapply(seq_len(nrow(prepared$historical_data)), function(i) {
      draws <- .extract_loo_draws(loo_posteriors[[i]], prepared$ids[i])
      sample_dat <- .historical_row_to_sample_dat(prepared$historical_data, i)
      metrics <- .loo_trial_metrics(draws, sample_dat)

      data.frame(
        trial_id = prepared$ids[i],
        observed_clinical = metrics$observed_clinical,
        posterior_mean_clinical = metrics$posterior_mean_clinical,
        squared_error = metrics$squared_error,
        predictive_percentile = metrics$predictive_percentile,
        predictive_lower_95 = metrics$predictive_lower_95,
        predictive_lower_90 = metrics$predictive_lower_90,
        predictive_upper_90 = metrics$predictive_upper_90,
        predictive_upper_95 = metrics$predictive_upper_95,
        covered_95 = metrics$covered_95,
        covered_90 = metrics$covered_90,
        stringsAsFactors = FALSE
      )
    })
    per_trial <- do.call(rbind, rows)
    rownames(per_trial) <- NULL

    summary <- data.frame(
      n_trials = nrow(per_trial),
      loo_cv_rmse = sqrt(mean(per_trial$squared_error)),
      coverage_95 = mean(per_trial$covered_95),
      coverage_90 = mean(per_trial$covered_90)
    )

    list(per_trial = per_trial, summary = summary)
  }

  .with_preserved_seed(seed, calculate())
}
