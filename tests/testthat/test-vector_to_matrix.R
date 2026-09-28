test_that("vector_to_matrix always returns one clinical and two surrogate components", {
  data("historical_posterior", package = "MBE")

  result <- MBE:::vector_to_matrix(historical_posterior[1:10, ])

  expect_length(result$mean, 10)
  expect_length(result$covar, 10)
  expect_equal(dim(result$mean[[1]]), c(3, 1))
  expect_equal(dim(result$covar[[1]]), c(3, 3))
})

test_that("vector_to_matrix rejects posterior draws for a third surrogate", {
  data("historical_posterior", package = "MBE")
  historical_posterior$muSur3 <- 0

  expect_error(
    MBE:::vector_to_matrix(historical_posterior[1:10, ]),
    "exactly two surrogates"
  )
})

test_that("vector_to_matrix supports fixed clinical intercept draws", {
  data("historical_posterior", package = "MBE")
  fixed_intercept_draws <- historical_posterior[1:10, ]
  fixed_intercept_draws$alphaCEonSur1Sur2 <- NULL

  result <- MBE:::vector_to_matrix(fixed_intercept_draws)

  expect_true(all(vapply(result$mean, function(x) is.finite(x[1]), logical(1))))
})
