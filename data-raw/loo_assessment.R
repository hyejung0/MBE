# Generate the compact package datasets from the completed leave-one-out run.
#
# This script requires the local, ignored directory containing all 66 posterior
# files and the saved assessment. The large posterior files remain in data-raw;
# only the two small derived tables are distributed with the package.

assessment_path <- file.path(
  "data-raw",
  "loo-random-inverse-gamma",
  "loo_assessment.rds"
)

if (!file.exists(assessment_path)) {
  stop(
    "The archived assessment is required; see the manuscript vignette for the current assessment workflow.",
    call. = FALSE
  )
}

loo_assessment <- readRDS(assessment_path)
loo_assessment_summary <- loo_assessment$summary
loo_assessment_by_trial <- loo_assessment$per_trial

# These tables preserve the original completed 66-fit assessment. Its coverage
# indicators used empirical percentiles and an unweighted surrogate update;
# the current prediction helper uses joint conditioning and surrogate-based
# weights, and coverage tests inclusive quantile bounds. Do not overwrite these archived results by
# treating them as a fresh calculation with the new interface.

stopifnot(
  nrow(loo_assessment_summary) == 1L,
  nrow(loo_assessment_by_trial) == 66L,
  loo_assessment_summary$n_trials == nrow(loo_assessment_by_trial),
  isTRUE(all.equal(
    loo_assessment_summary$loo_cv_rmse,
    sqrt(mean(loo_assessment_by_trial$squared_error))
  )),
  isTRUE(all.equal(
    loo_assessment_summary$coverage_95,
    mean(loo_assessment_by_trial$covered_95)
  )),
  isTRUE(all.equal(
    loo_assessment_summary$coverage_90,
    mean(loo_assessment_by_trial$covered_90)
  ))
)

usethis::use_data(
  loo_assessment_summary,
  loo_assessment_by_trial,
  overwrite = TRUE,
  compress = "xz"
)
