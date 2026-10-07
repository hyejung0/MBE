#' Calculate RMSE from Paired Observed and Estimated Clinical Effects
#'
#' @description Calculates root mean squared error from user-supplied pairs.
#' Model fitting and selection of the point estimate are separate steps. No
#' historical models are fitted and no MBE updates are performed here.
#'
#' @param data A numeric matrix or data frame with one row per trial and exactly
#'   two columns named `observed` and `estimated`, in either order. `observed`
#'   contains the observed clinical treatment-effect estimates. `estimated`
#'   contains user-selected point estimates of the corresponding true clinical
#'   effects, for example `MBE(...)$post_mean[1, 1]`. Values must be finite;
#'   missing values are not silently omitted. Pair trials by ID before assembling
#'   this table; the function uses the pairs within each row.
#'
#' @details The result is `sqrt(mean((observed - estimated)^2))`, in the same
#' units as the input effects. It is a LOO-CV measure only when the estimates
#' were obtained using a historical model fitted without the corresponding
#' trial. The function cannot verify how the estimates were obtained.
#'
#' An MBE posterior mean estimates a true effect; it is not the known true
#' effect. When an MBE update incorporates the held-out clinical estimate,
#' this RMSE measures agreement with that observed estimate after updating.
#' It does not measure prediction from surrogates alone or error against a
#' known simulation truth. The intercept, priors, and MBE update settings are
#' chosen by the user before this calculation.
#'
#' @return One numeric RMSE value.
#' @export
#' @examples
#' pairs <- data.frame(observed = c(-0.2, 0.1), estimated = c(-0.1, 0.05))
#' loo_cv_rmse(pairs)
#' loo_cv_rmse(as.matrix(pairs))
#'
#' data("loo_assessment_by_trial")
#' loo_cv_rmse(data.frame(
#'   observed = loo_assessment_by_trial$observed_clinical,
#'   estimated = loo_assessment_by_trial$posterior_mean_clinical
#' ))
loo_cv_rmse <- function(data) {
  if ((!is.data.frame(data) && !is.matrix(data)) ||
      nrow(data) < 1L || ncol(data) != 2L ||
      !setequal(colnames(data), c("observed", "estimated"))) {
    stop(
      "`data` must be a non-empty two-column matrix or data frame named `observed` and `estimated`.",
      call. = FALSE
    )
  }
  pairs <- as.data.frame(data)
  valid <- vapply(pairs, function(x) {
    is.numeric(x) && is.null(dim(x)) && all(is.finite(x))
  }, logical(1))
  if (!all(valid)) {
    stop("Both columns must contain only finite numeric values.", call. = FALSE)
  }

  # Scaling avoids overflow when finite errors are large enough that squaring
  # them directly would exceed the floating-point range.
  error <- pairs$observed - pairs$estimated
  scale <- max(abs(error))
  if (scale == 0 || is.infinite(scale)) {
    return(scale)
  }
  scale * sqrt(mean((error / scale)^2))
}
