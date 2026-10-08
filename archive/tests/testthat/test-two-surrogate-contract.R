new_trial_input <- function(trial_sim_dat) {
  list(
    ClnEst = trial_sim_dat$CE_est[1],
    ClnSE = trial_sim_dat$CE_se[1],
    Sur1Est = trial_sim_dat$Sur1_est[1],
    Sur1SE = trial_sim_dat$Sur1_se[1],
    Sur2Est = trial_sim_dat$Sur2_est[1],
    Sur2SE = trial_sim_dat$Sur2_se[1],
    R1Clin = trial_sim_dat$Cor_CE_Sur1[1],
    R2Clin = trial_sim_dat$Cor_CE_Sur2[1],
    R12 = trial_sim_dat$Cor_Sur1_Sur2[1]
  )
}

historical_input <- function(trial_sim_dat) {
  as.data.frame(trial_sim_dat)[, c(
    "CE_est", "CE_se",
    "Sur1_est", "Sur1_se",
    "Sur2_est", "Sur2_se",
    "Cor_CE_Sur1", "Cor_CE_Sur2", "Cor_Sur1_Sur2"
  ), drop = FALSE]
}

test_that("MBE returns three components in the documented order", {
  data("historical_posterior", package = "MBE")
  data("trial_sim_dat", package = "MBE")

  set.seed(1)
  result <- MBE(
    historical_posterior[1:25, ],
    new_trial_input(trial_sim_dat)
  )

  expect_equal(rownames(result$post_mean), c("clinical", "surrogate1", "surrogate2"))
  expect_equal(dim(result$post_var), c(3, 3))
  expect_equal(colnames(result$post_psi0), c("clinical", "surrogate1", "surrogate2"))
  expect_named(
    result,
    c(
      "post_mean", "post_var",
      "post_quantiles_clinical",
      "post_quantiles_surrogate1",
      "post_quantiles_surrogate2",
      "post_psi0", "weight_data", "importance_diagnostics"
    )
  )
  expect_named(
    result$importance_diagnostics,
    c(
      "effective_sample_size",
      "relative_effective_sample_size",
      "maximum_normalized_weight"
    )
  )
  expect_gte(result$importance_diagnostics[["effective_sample_size"]], 1)
  expect_lte(result$importance_diagnostics[["effective_sample_size"]], 25)
  expect_gt(result$importance_diagnostics[["relative_effective_sample_size"]], 0)
  expect_lte(result$importance_diagnostics[["relative_effective_sample_size"]], 1)
})

test_that("MBE rejects a third surrogate endpoint", {
  data("historical_posterior", package = "MBE")
  data("trial_sim_dat", package = "MBE")

  sample_dat <- new_trial_input(trial_sim_dat)
  sample_dat$Sur3Est <- 0

  expect_error(
    MBE(historical_posterior[1:10, ], sample_dat),
    "exactly the nine fields"
  )
})

test_that("historical fitting input is restricted to two surrogates", {
  data("trial_sim_dat", package = "MBE")
  data <- historical_input(trial_sim_dat)
  data$Sur3_est <- 0

  expect_error(
    historical_model_fit_2surrogates(data, nchains = 1, ncores = 1, niter = 1, nwarmup = 1),
    "exactly the nine named columns"
  )
})

test_that("historical fit summaries pass measures individually", {
  draws <- posterior::as_draws_df(data.frame(theta = seq_len(10)))
  fake_fit <- list(
    summary = function(variables, ...) {
      posterior::summarise_draws(draws, ...)
    }
  )

  result <- MBE:::.summarize_historical_fit(fake_fit, "theta")

  expect_named(
    result,
    c("variable", "mean", "median", "sd", "mad", "q2.5", "q97.5")
  )
})

test_that("packaged datasets contain no third-surrogate fields", {
  data("trial_sim_dat", package = "MBE")

  expect_false(any(grepl("Sur3", names(trial_sim_dat))))
})
