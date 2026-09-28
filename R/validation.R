.mbe_sample_fields <- c(
  "ClnEst", "ClnSE",
  "Sur1Est", "Sur1SE",
  "Sur2Est", "Sur2SE",
  "R1Clin", "R2Clin", "R12"
)

.mbe_historical_fields <- c(
  "CE_est", "CE_se",
  "Sur1_est", "Sur1_se",
  "Sur2_est", "Sur2_se",
  "Cor_CE_Sur1", "Cor_CE_Sur2", "Cor_Sur1_Sur2"
)

.mbe_mcmc_parameters <- c(
  "b1CEonSur1Sur2",
  "b2CEonSur1Sur2",
  "SigSqCEonSur1Sur2",
  "alphaSur1onSur2",
  "bSur1onSur2",
  "SigSqSur1onSur2",
  "muSur2",
  "sigSqSur2"
)

.format_field_list <- function(x) {
  paste(x, collapse = ", ")
}

.validate_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
}

.validate_positive_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be one finite number greater than zero.", name), call. = FALSE)
  }
}

.validate_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != as.integer(x)) {
    stop(sprintf("`%s` must be one positive integer.", name), call. = FALSE)
  }
}

.validate_correlation <- function(x, name) {
  if (any(x < -1 | x > 1)) {
    stop(sprintf("`%s` must contain values between -1 and 1.", name), call. = FALSE)
  }
}

.validate_positive_definite <- function(x, label) {
  is_positive_definite <- tryCatch({
    chol(x)
    TRUE
  }, error = function(e) FALSE)

  if (!is_positive_definite) {
    stop(sprintf("%s must be positive definite.", label), call. = FALSE)
  }
}

.validate_sample_dat <- function(sample_dat) {
  if (!is.list(sample_dat) || is.null(names(sample_dat))) {
    stop("`sample_dat` must be a named list for one clinical and exactly two surrogate endpoints.", call. = FALSE)
  }

  missing_fields <- setdiff(.mbe_sample_fields, names(sample_dat))
  extra_fields <- setdiff(names(sample_dat), .mbe_sample_fields)
  if (length(missing_fields) || length(extra_fields)) {
    details <- character()
    if (length(missing_fields)) {
      details <- c(details, paste0("missing: ", .format_field_list(missing_fields)))
    }
    if (length(extra_fields)) {
      details <- c(details, paste0("unexpected: ", .format_field_list(extra_fields)))
    }
    stop(
      paste0(
        "`sample_dat` must contain exactly the nine fields for one clinical and two surrogate endpoints (",
        paste(details, collapse = "; "), ")."
      ),
      call. = FALSE
    )
  }

  scalar_numeric <- vapply(
    sample_dat[.mbe_sample_fields],
    function(x) is.numeric(x) && length(x) == 1L && is.finite(x),
    logical(1)
  )
  if (!all(scalar_numeric)) {
    stop(
      paste0(
        "Every `sample_dat` field must be one finite numeric value. Invalid fields: ",
        .format_field_list(names(scalar_numeric)[!scalar_numeric]), "."
      ),
      call. = FALSE
    )
  }

  if (sample_dat$ClnSE <= 0 || sample_dat$Sur1SE <= 0 || sample_dat$Sur2SE <= 0) {
    stop("`ClnSE`, `Sur1SE`, and `Sur2SE` must be greater than zero.", call. = FALSE)
  }

  for (field in c("R1Clin", "R2Clin", "R12")) {
    .validate_correlation(sample_dat[[field]], field)
  }

  invisible(sample_dat)
}

.validate_mcmc_dat <- function(mcmc_dat) {
  if (!is.data.frame(mcmc_dat) && !is.matrix(mcmc_dat)) {
    stop("`mcmc_dat` must be a data frame or matrix of posterior draws.", call. = FALSE)
  }
  if (!nrow(mcmc_dat)) {
    stop("`mcmc_dat` must contain at least one posterior draw.", call. = FALSE)
  }

  missing_parameters <- setdiff(.mbe_mcmc_parameters, colnames(mcmc_dat))
  if (length(missing_parameters)) {
    stop(
      paste0(
        "`mcmc_dat` is missing parameters required by the two-surrogate model: ",
        .format_field_list(missing_parameters), "."
      ),
      call. = FALSE
    )
  }

  third_surrogate_parameters <- grep("sur3|surrogate3", colnames(mcmc_dat), value = TRUE, ignore.case = TRUE)
  if (length(third_surrogate_parameters)) {
    stop(
      paste0(
        "`mcmc_dat` contains third-surrogate parameters, but MBE supports exactly two surrogates: ",
        .format_field_list(third_surrogate_parameters), "."
      ),
      call. = FALSE
    )
  }

  columns_to_check <- intersect(
    c("alphaCEonSur1Sur2", .mbe_mcmc_parameters),
    colnames(mcmc_dat)
  )
  required_values <- as.matrix(as.data.frame(mcmc_dat)[, columns_to_check, drop = FALSE])
  if (!is.numeric(required_values) || any(!is.finite(required_values))) {
    stop("The required columns in `mcmc_dat` must contain only finite numeric values.", call. = FALSE)
  }
  if (any(required_values[, c("SigSqCEonSur1Sur2", "SigSqSur1onSur2", "sigSqSur2"), drop = FALSE] <= 0)) {
    stop("Variance parameters in `mcmc_dat` must be greater than zero.", call. = FALSE)
  }

  invisible(mcmc_dat)
}

.validate_historical_data <- function(data) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data frame or matrix.", call. = FALSE)
  }
  if (!nrow(data)) {
    stop("`data` must contain at least one trial.", call. = FALSE)
  }
  if (is.null(colnames(data))) {
    stop("`data` must have named columns for one clinical and exactly two surrogate endpoints.", call. = FALSE)
  }

  missing_fields <- setdiff(.mbe_historical_fields, colnames(data))
  extra_fields <- setdiff(colnames(data), .mbe_historical_fields)
  if (length(missing_fields) || length(extra_fields)) {
    details <- character()
    if (length(missing_fields)) {
      details <- c(details, paste0("missing: ", .format_field_list(missing_fields)))
    }
    if (length(extra_fields)) {
      details <- c(details, paste0("unexpected: ", .format_field_list(extra_fields)))
    }
    stop(
      paste0(
        "`data` must contain exactly the nine named columns for one clinical and two surrogate endpoints (",
        paste(details, collapse = "; "), ")."
      ),
      call. = FALSE
    )
  }

  values <- as.matrix(as.data.frame(data)[, .mbe_historical_fields, drop = FALSE])
  if (!is.numeric(values) || any(!is.finite(values))) {
    stop("All columns in `data` must contain only finite numeric values.", call. = FALSE)
  }
  if (any(values[, c("CE_se", "Sur1_se", "Sur2_se"), drop = FALSE] <= 0)) {
    stop("`CE_se`, `Sur1_se`, and `Sur2_se` must be greater than zero.", call. = FALSE)
  }
  for (field in c("Cor_CE_Sur1", "Cor_CE_Sur2", "Cor_Sur1_Sur2")) {
    .validate_correlation(values[, field], field)
  }

  invisible(data)
}
