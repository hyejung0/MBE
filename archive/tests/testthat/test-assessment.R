prediction_test_input <- function() {
  list(ClnSE = 1, Sur1Est = 0.2, Sur1SE = 1, Sur2Est = -0.1, Sur2SE = 1,
       R1Clin = 0.2, R2Clin = 0.1, R12 = 0.1)
}

test_that("RMSE accepts paired tables and matrices with columns in either order", {
  pairs <- data.frame(observed = c(1, 4), estimated = c(2, 2))
  expect_equal(loo_cv_rmse(pairs), sqrt(2.5))
  expect_equal(loo_cv_rmse(as.matrix(pairs[, 2:1])), sqrt(2.5))
  expect_equal(loo_cv_rmse(data.table::as.data.table(pairs)), sqrt(2.5))
  expect_equal(loo_cv_rmse(data.frame(observed = 2, estimated = 2)), 0)
  expect_equal(loo_cv_rmse(data.frame(observed = 3, estimated = 1)), 2)
  expect_equal(loo_cv_rmse(data.frame(observed = 1e200, estimated = 0)), 1e200)
})

test_that("RMSE rejects ambiguous pairs and invalid values", {
  for (x in list(matrix(1:4, 2), data.frame(observed = 1),
                 data.frame(observed = numeric(), estimated = numeric()),
                 data.frame(observed = 1, estimated = 2, trial_id = 1))) {
    expect_error(loo_cv_rmse(x), "two-column")
  }
  for (bad in list(NA_real_, Inf, "1", factor("1"))) {
    expect_error(loo_cv_rmse(data.frame(observed = 1, estimated = bad)), "finite numeric")
  }
})

test_that("coverage calculates user-selected intervals in the requested order", {
  calls <- 0L
  testthat::local_mocked_bindings(
    clinical_predictive_distribution = function(...) { calls <<- calls + 1L; 0:100 },
    .package = "MBE"
  )
  result <- tolerance_interval_coverage(94, NULL, NULL,
    rbind(c(.075, .925), c(.025, .975), c(.05, .95)))
  expect_equal(result$level, c(85, 95, 90))
  expect_equal(result$lower, c(7.5, 2.5, 5))
  expect_equal(result$upper, c(92.5, 97.5, 95))
  expect_identical(result$covered, c(FALSE, TRUE, TRUE))
  expect_equal(calls, 1L)
  expect_equal(nrow(tolerance_interval_coverage(50, NULL, NULL)), 1L)
  expect_equal(tolerance_interval_coverage(50, NULL, NULL)$level, 95)
  asymmetric <- tolerance_interval_coverage(50, NULL, NULL, c(0, .9))
  expect_equal(asymmetric$lower, 0)
  expect_equal(asymmetric$upper, 90)
})

test_that("coverage includes interval boundaries and handles ties", {
  bounds <- matrix(c(.05, .95), nrow = 1)
  capture <- MBE:::.predictive_interval_capture
  expect_true(capture(5, 0:100, bounds)$covered)
  expect_true(capture(95, 0:100, bounds)$covered)
  expect_false(capture(4.99, 0:100, bounds)$covered)
  expect_true(all(capture(2, rep(2, 10), bounds)$covered))
  expect_false(any(capture(3, rep(2, 10), bounds)$covered))
})

test_that("coverage validates observations and quantile probability pairs", {
  for (x in list(NA_real_, Inf, c(1, 2), "1", numeric())) {
    expect_error(tolerance_interval_coverage(x, NULL, NULL), "observed")
  }
  for (x in list(numeric(), NA, Inf, 95, c(.9, .1), c(.5, .5), c(-.1, .9),
                 c(.1, 1.1), c(.025, NA), c(".025", ".975"), matrix(1:6, 2))) {
    expect_error(tolerance_interval_coverage(1, NULL, NULL, x), "probs")
  }
})

test_that("assessment functions do not fit models or call MBE", {
  testthat::local_mocked_bindings(
    historical_model_fit_2surrogates = function(...) stop("Unexpected refitting"),
    MBE = function(...) stop("Unexpected MBE update"),
    .package = "MBE"
  )
  set.seed(42)
  before <- .Random.seed
  expect_equal(loo_cv_rmse(data.frame(observed = 1, estimated = 2)), 1)
  expect_identical(.Random.seed, before)
  data("historical_posterior", package = "MBE")
  sampling <- prediction_test_input()
  result <- tolerance_interval_coverage(0, historical_posterior[1:20, ],
                                        sampling, seed = 42)
  expect_type(result$covered, "logical")
  expect_identical(.Random.seed, before)
})

test_that("predictive components match analytic normal conditioning", {
  prior <- list(mean = list(c(1, 2, 3)), covar = list(diag(3)))
  sampling <- matrix(c(4, 1, 0, 1, 1, 0, 0, 0, 1), 3)
  result <- MBE:::.clinical_predictive_components(prior, sampling, c(4, 3))
  expect_equal(result$mean, 2)
  expect_equal(result$variance, 4.5)
  expect_equal(result$weight, 1)
  # With independent sampling errors and independent model components, clinical
  # predictive variance is residual variance + clinical sampling variance.
  independent <- MBE:::.clinical_predictive_components(prior, diag(c(4, 1, 1)), c(4, 3))
  expect_equal(independent$mean, 1)
  expect_equal(independent$variance, 5)

  prior$mean[[2]] <- c(1, 20, 3)
  prior$covar[[2]] <- diag(3)
  weighted <- MBE:::.clinical_predictive_components(prior, sampling, c(4, 3))
  expect_gt(weighted$weight[1], .999)
  expect_equal(sum(weighted$weight), 1)
})

test_that("prediction is reproducible and ignores held-out clinical observations", {
  data("historical_posterior", package = "MBE")
  draws <- historical_posterior[1:30, ]
  input <- prediction_test_input()
  set.seed(234)
  before <- .Random.seed
  predicted <- clinical_predictive_distribution(draws, input, n_draws = 200, seed = 42)
  expect_length(predicted, 200)
  expect_true(all(is.finite(predicted)))
  expect_identical(.Random.seed, before)
  input$ClnEst <- 1e6
  expect_identical(predicted, clinical_predictive_distribution(draws, input, 200, 42))
  input$ClnEst <- NA_real_
  expect_identical(predicted, clinical_predictive_distribution(draws, input, 200, 42))
  actual <- tolerance_interval_coverage(0, draws, input, n_draws = 200, seed = 42)
  expected <- stats::quantile(predicted, c(.025, .975), names = FALSE)
  expect_equal(c(actual$lower, actual$upper), expected)
})

test_that("prediction supports fixed intercepts and user-selected carryover", {
  data("historical_posterior", package = "MBE")
  draws <- as.data.frame(historical_posterior[1:20, ])
  draws$alphaCEonSur1Sur2 <- NULL
  fixed <- clinical_predictive_distribution(draws, prediction_test_input(), seed = 42)
  draws$alphaCEonSur1Sur2 <- 0
  expect_identical(fixed, clinical_predictive_distribution(draws, prediction_test_input(), seed = 42))
  expect_true(all(is.finite(clinical_predictive_distribution(
    draws, prediction_test_input(), diffuse = TRUE, intercept0 = TRUE, seed = 42
  ))))
})

test_that("prediction rejects invalid inputs before simulation", {
  data("historical_posterior", package = "MBE")
  draws <- historical_posterior[1:10, ]
  input <- prediction_test_input()
  for (n in list(0, NA, Inf, 1.5, 1e20)) {
    expect_error(clinical_predictive_distribution(draws, input, n_draws = n), "n_draws")
    expect_error(clinical_predictive_distribution(draws, input, seed = n), "seed")
  }
  input$ClnSE <- -1
  expect_error(clinical_predictive_distribution(draws, input), "greater than zero")
  input <- prediction_test_input()
  input$R1Clin <- .99
  input$R2Clin <- -.99
  expect_error(clinical_predictive_distribution(draws, input), "positive definite")
})

test_that("archived LOO tables retain their values under the new RMSE interface", {
  data("loo_assessment_summary", package = "MBE")
  data("loo_assessment_by_trial", package = "MBE")
  expect_equal(nrow(loo_assessment_by_trial), loo_assessment_summary$n_trials)
  expect_equal(loo_cv_rmse(data.frame(
    observed = loo_assessment_by_trial$observed_clinical,
    estimated = loo_assessment_by_trial$posterior_mean_clinical
  )), loo_assessment_summary$loo_cv_rmse)
  expect_equal(mean(loo_assessment_by_trial$covered_95), loo_assessment_summary$coverage_95)
  expect_equal(mean(loo_assessment_by_trial$covered_90), loo_assessment_summary$coverage_90)
})
