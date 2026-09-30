make_loo_draws <- function(n = 100L) {
  data.frame(
    alphaCEonSur1Sur2 = rep(0, n),
    b1CEonSur1Sur2 = rep(0.4, n),
    b2CEonSur1Sur2 = rep(0.2, n),
    SigSqCEonSur1Sur2 = rep(0.25, n),
    alphaSur1onSur2 = rep(0, n),
    bSur1onSur2 = rep(0.1, n),
    SigSqSur1onSur2 = rep(0.5, n),
    muSur2 = rep(0, n),
    sigSqSur2 = rep(0.5, n)
  )
}

test_that("clinical predictive SD includes residual and sampling variance", {
  expect_equal(
    MBE:::.clinical_predictive_sd(c(4, 9), 3),
    sqrt(c(13, 18))
  )
})

test_that("LOO assessment returns per-trial and aggregate manuscript metrics", {
  data("trial_sim_dat", package = "MBE")
  held_out <- as.data.frame(trial_sim_dat[1:2, ])
  draws <- list(`1` = make_loo_draws(), `2` = make_loo_draws())

  result <- loo_cv_model_assessment(draws, held_out, seed = 42)
  repeated <- loo_cv_model_assessment(draws, held_out, seed = 42)

  expect_equal(result, repeated)
  expect_equal(result$per_trial$trial_id, c("1", "2"))
  expect_type(result$per_trial$covered_95, "logical")
  expect_type(result$per_trial$covered_90, "logical")
  expect_equal(
    result$summary$loo_cv_rmse,
    sqrt(mean(result$per_trial$squared_error))
  )
  expect_equal(result$summary$coverage_95, mean(result$per_trial$covered_95))
  expect_equal(result$summary$coverage_90, mean(result$per_trial$covered_90))
})

test_that("LOO assessment matches named posteriors to trial IDs", {
  data("trial_sim_dat", package = "MBE")
  held_out <- as.data.frame(trial_sim_dat[1:2, ])
  draws <- list(`2` = make_loo_draws(), `1` = make_loo_draws())

  result <- loo_cv_model_assessment(draws, held_out, seed = 42)

  expect_equal(result$per_trial$trial_id, c("1", "2"))
})

test_that("LOO assessment rejects mismatched posterior inputs", {
  data("trial_sim_dat", package = "MBE")
  held_out <- as.data.frame(trial_sim_dat[1:2, ])

  expect_error(
    loo_cv_model_assessment(list(make_loo_draws()), held_out),
    "1 elements.*2 trials"
  )
  expect_error(
    loo_cv_model_assessment(
      list(first = make_loo_draws(), second = make_loo_draws()),
      held_out
    ),
    "matching every trial ID"
  )
})

test_that("LOO fit runner saves compact restartable records", {
  data("trial_sim_dat", package = "MBE")
  historical_data <- as.data.frame(trial_sim_dat[1:2, ])
  mock_draws <- make_loo_draws(10)
  fit_count <- 0L
  testthat::local_mocked_bindings(
    historical_model_fit_2surrogates = function(...) {
      fit_count <<- fit_count + 1L
      list(posterior_samples = mock_draws)
    },
    .package = "MBE"
  )

  output_dir <- tempfile("mbe-loo-fits-")
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  paths <- fit_loo_historical_models(
    historical_data,
    output_dir = output_dir,
    seed = 10,
    show_messages = FALSE
  )

  expect_equal(fit_count, 2L)
  expect_true(all(file.exists(paths)))
  first_record <- readRDS(paths[[1]])
  expect_equal(first_record$held_out_trial, "1")
  expect_equal(first_record$specification$random_intercept, TRUE)
  expect_equal(first_record$specification$prior_for_uncertainty, "inverse_gamma")
  expect_named(first_record$posterior_draws, MBE:::.mbe_loo_draw_fields)

  repeated_paths <- fit_loo_historical_models(
    historical_data,
    output_dir = output_dir,
    seed = 10,
    show_messages = FALSE
  )
  expect_equal(repeated_paths, paths)
  expect_equal(fit_count, 2L)

  expect_error(
    fit_loo_historical_models(
      historical_data,
      output_dir = output_dir,
      seed = 11,
      show_messages = FALSE
    ),
    "does not match the requested trial, seed, or model specification"
  )

  assessment <- loo_cv_model_assessment(output_dir, historical_data, seed = 42)
  expect_equal(assessment$per_trial$trial_id, c("1", "2"))
})

test_that("packaged LOO assessment tables are internally consistent", {
  data("loo_assessment_summary", package = "MBE")
  data("loo_assessment_by_trial", package = "MBE")

  expect_equal(nrow(loo_assessment_summary), 1L)
  expect_equal(nrow(loo_assessment_by_trial), 66L)
  expect_equal(
    loo_assessment_summary$n_trials,
    nrow(loo_assessment_by_trial)
  )
  expect_equal(
    loo_assessment_summary$loo_cv_rmse,
    sqrt(mean(loo_assessment_by_trial$squared_error))
  )
  expect_equal(
    loo_assessment_summary$coverage_95,
    mean(loo_assessment_by_trial$covered_95)
  )
  expect_equal(
    loo_assessment_summary$coverage_90,
    mean(loo_assessment_by_trial$covered_90)
  )
})
