# Internal session cache for CmdStanModel objects
.pkg_models <- new.env(parent = emptyenv())

#' Get or Compile a CmdStan Model
#'
#' @param model_name Name of the Stan file (without .stan extension).
#' @return A `CmdStanModel` object.
#' @noRd
get_cmdstan_model <- function(model_name) {
  # 1. Return in-memory model if already loaded in this session
  if (exists(model_name, envir = .pkg_models)) {
    return(get(model_name, envir = .pkg_models))
  }

  # 2. Locate the .stan file in the installed package
  stan_file <- system.file("stan", paste0(model_name, ".stan"), package = "MBE")
  if (stan_file == "" || !file.exists(stan_file)) {
    stop(sprintf("Stan file '%s.stan' not found in package.", model_name), call. = FALSE)
  }

  # 3. Set up a persistent cache folder for the compiled C++ binary
  cache_dir <- tools::R_user_dir("MBE", which = "cache")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  # 4. Compile or load existing executable from cache
  mod <- cmdstanr::cmdstan_model(
    stan_file = stan_file,
    dir = cache_dir
  )

  # 5. Cache in memory for instant reuse
  assign(model_name, mod, envir = .pkg_models)
  return(mod)
}

#' Trial-level Meta-Analysis on Historical Data using 2-stage Random Effects Bayesian Model (2 surrogates)
#'
#' @description
#' Fits a trial-level meta-analysis model using a 2-stage random effects Bayesian model for two surrogate endpoints and a definitive clinical endpoint, as described in Lee et al. (2026).
#'
#' @param data n by p data set with n trials and p covariates.
#'   The first column should be the estimated treatment effect on the definitive
#'   clinical endpoint. The following column should be standard error of the treatment effect on the clinical endpoint.
#'   The next two columns should be the estimated treatment effect on the first surrogate endpoint and its standard error.
#'   The next two columns should be the estimated treatment effect on the second surrogate endpoint and its standard error.
#'   The last three columns should be the correlation between the clinical endpoint
#'   and the first surrogates, in the order of clinical endpoint and first surrogate,
#'   clinical endpoint and second surrogate, and first surrogate and second surrogate.
#' @param random_intercept logical indicating whether to include intercept in
#'   regression modeling.
#' @param prior_for_uncertainty character indicating the prior distribution for the uncertainty parameters.
#'   Options are "half_normal" for half-normal prior on the standard deviation
#'   scale and "inverse_gamma" for inverse-gamma prior on the variance scale.
#'   Default is "half_normal".
#' @param nchains number of chains for MCMC sampling. Default to 4.
#' @param ncores number of chains to run in parallel. Default to 1
#'   (sequential processing). If set to a value greater than 1,
#'   the function will use parallel processing to run multiple chains simultaneously.
#' @param niter number of post-warmup samples to save for each chain. Default is 2000.
#' @param nwarmup number of warmup iterations for MCMC sampling.
#' @param output_dir Character string. The folder path where CmdStan should save
#'   the raw MCMC chain CSV files. Defaults to \code{tempdir()}, which saves files
#'   to a temporary session folder that is automatically deleted when R closes.
#'   To keep the CSV files permanently, provide a local directory path (e.g., \code{"./"}
#'   for the current working directory, or \code{"./mcmc_output"}).
#' @param show_messages Logical; whether to print progress to console. Default is TRUE.
#' @param ... Additional arguments forwarded to \code{cmdstanr}'s \code{$sample()}
#'   method (e.g., \code{adapt_delta}, \code{max_treedepth}, \code{seed}, \code{refresh}).
#'   The \code{cmdstanr}'s default values are used for \code{adapt_delta = 0.9} and \code{max_treedepth = 12}.
#'   In addition, \code{chains}, \code{parallel_chains}, \code{iter_sampling}, and \code{iter_warmup} are set to the values of \code{nchains}, \code{ncores}, \code{niter}, and \code{nwarmup}, respectively, unless overridden in \code{...}.
#'
#' @details
#' This function fits a trial-level meta-analysis model using a 2-stage random effects Bayesian model
#' under chronic kidney disease (CKD) context first introduced by Lee et al. (2026).
#' For a true treatment effect on clinical endpoint \eqn{\theta_i}, and true treatment effect on two surrogate endpoints
#' \eqn{\gamma_{i, 1}} and \eqn{\gamma_{i, 2}} for the \eqn{i}-th trial, the model is specified as follows:
#' \deqn{\hat{\psi}_i \mid \psi_i \sim N_3(\psi_i, \Sigma_{y, i}), \quad \psi_i = (\mu, \Sigma)^T}{hat(psi)_i | psi_i ~ N_3(psi_i, Sigma_{y,i}), psi_i = (\mu, \Sigma)^T}
#' where \eqn{\psi_i = (\theta_i, \gamma_{i, 1}, \gamma_{i, 2})^T}{\psi_i = (\theta_i, \gamma_{i,1}, \gamma_{i,2})^T} is the true
#' treatment effect vector for the \eqn{i}-th trial, \eqn{\hat{\psi}_i}{hat(\psi)_i} is
#' the observed (estimated) treatment effect vector for the \eqn{i}-th trial,
#' and \eqn{\Sigma_{y, i}}{\Sigma_{y,i}} is the covariance matrix of the corresponding \eqn{\hat{\psi}_i}{hat(\psi)_i}.
#'
#' The density of \eqn{\psi}{psi} can be expressed as a sequence of conditional densities:
#' \deqn{\gamma_{i, 2} \sim N(\mu_2, \sigma_2^2)}{\gamma_{i,2} ~ N(\mu_2, \sigma_2^2)}
#' \deqn{\gamma_{i, 1} \mid \gamma_{i, 2} \sim N(\alpha_{\gamma_1} + \beta_{\gamma_1} \cdot \gamma_{i, 2}, \lambda^2_{\gamma_1})}{\gamma_{i,1} | \gamma_{i,2} ~ N(\alpha_{\gamma 1} + \beta_{\gamma 1} * \gamma_{i,2}, \lambda^2_{\gamma 1})}
#' \deqn{\theta_i \mid \gamma_{i, 1}, \gamma_{i, 2} \sim N(\alpha_\theta + \beta_{\gamma_1} \cdot \gamma_{i, 1} + \beta_{\gamma_2} \cdot \gamma_{i, 2}, \lambda^2_\theta)}{\theta_i | \gamma_{i,1}, \gamma_{i,2} ~ N(\alpha_\theta + \beta_{\gamma 1} * \gamma_{i,1} + \beta_{\gamma 2} * \gamma_{i,2}, \lambda^2_\theta)}
#'
#' This function fits the model using MCMC sampling on historical RCTs and returns
#' the posterior distribution of the model parameters.
#'
#' @return A list containing the fitted model object and other relevant information:
#' \describe{
#'   \item{diagnostics}{Diagnostic summary table from \code{fit$diagnostic_summary()}.}
#'   \item{posterior_samples}{Sampled posterior distribution draws as a data frame.}
#'   \item{summary}{A summary of the posterior distribution of parameters of interest.}
#'   \item{loo}{Efficient approximate leave-one-out (LOO) cross-validation from \code{loo::loo()}.}
#'   \item{waic}{Widely applicable information criterion (WAIC) calculated with \code{loo::waic()}.}
#'   \item{fit}{The fitted \code{CmdStanMCMC} object from \code{cmdstanr}.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' data("trial_sim_dat", package = "MBE")
#'
#' fit <- historical_model_fit_2surrogates(
#'   data = trial_sim_dat,
#'   random_intercept = TRUE,
#'   nchains = 4,
#'   ncores = 4,
#'   niter = 2000,
#'   nwarmup = 1000,
#'   adapt_delta = 0.95,  # Target acceptance rate (resolves divergences)
#'   max_treedepth = 15,  # Tree depth limit (resolves treedepth saturation)
#'   refresh = 250        # Print updates every 250 iterations
#' )
#' }
historical_model_fit_2surrogates<-function(data,
                                           random_intercept=TRUE,
                                           prior_for_uncertainty="half_normal",
                                           output_dir = tempdir(),
                                           show_messages=TRUE,
                                           nchains=4,
                                           ncores=1,
                                           niter=2000,
                                           nwarmup=1000,
                                           ...){

  # Check data format
  if (is.null(data) || ncol(data) != 9) {
    stop("Data must have exactly 9 columns. Please see help file for details.", call. = FALSE)
  }

  # Choose Stan file to load
  if (isTRUE(random_intercept) && prior_for_uncertainty == "half_normal") {
    model_name <- "random_beta0_halfNormal"
  } else if (isTRUE(random_intercept) && prior_for_uncertainty == "inverse_gamma") {
    model_name <- "random_beta0_invGamma"
  } else if (!isTRUE(random_intercept) && prior_for_uncertainty == "half_normal") {
    model_name <- "fix_beta0_halfNormal"
  } else if (!isTRUE(random_intercept) && prior_for_uncertainty == "inverse_gamma") {
    model_name <- "fix_beta0_invGamma"
  } else {
    stop(
      "Invalid combination of random_intercept and prior_for_uncertainty. ",
      "random_intercept must be TRUE or FALSE, and prior_for_uncertainty must be 'half_normal' or 'inverse_gamma'.",
      call. = FALSE
    )
  }

  #Make sure it is data.table
  dt <- data.table::as.data.table(data.table::copy(data))
  #standardize column names on a local copy
  data.table::setnames(dt, c("CE_est", "CE_se", "Sur1_est", "Sur1_se", "Sur2_est", "Sur2_se", "Cor_CE_Sur1", "Cor_CE_Sur2", "Cor_Sur1_Sur2"))

  N <- nrow(dt)
  obs_mean <- vector("list", N)
  obs_var  <- vector("list", N)

  # Construct Stan data structures (pure numeric vectors and base matrices)
  for (j in seq_len(N)) {
    obs_mean[[j]] <- c(dt$CE_est[j], dt$Sur1_est[j], dt$Sur2_est[j])

    v1  <- dt$CE_se[j]^2
    v2  <- dt$Sur1_se[j]^2
    v3  <- dt$Sur2_se[j]^2
    c12 <- dt$Cor_CE_Sur1[j] * dt$CE_se[j] * dt$Sur1_se[j]
    c13 <- dt$Cor_CE_Sur2[j] * dt$CE_se[j] * dt$Sur2_se[j]
    c23 <- dt$Cor_Sur1_Sur2[j] * dt$Sur1_se[j] * dt$Sur2_se[j]

    obs_var[[j]] <- matrix(
      c(
        v1,  c12, c13,
        c12, v2,  c23,
        c13, c23, v3
      ),
      nrow = 3, ncol = 3, byrow = TRUE
    )
  }


  #collect them together as a more comprehensive list
  data_sub = list(
    N=N,
    obs_mean = obs_mean,
    obs_var = obs_var
  )


  # Retrieve compiled model
  mod <- get_cmdstan_model(model_name)


  #Set MCMC parameters
  default_args <- list(
    data = data_sub,
    output_dir = output_dir,
    show_messages = show_messages,
    show_exceptions = TRUE,
    chains = nchains,
    parallel_chains = ncores,
    iter_sampling = niter,
    iter_warmup = nwarmup,
    adapt_delta = 0.9,
    max_treedepth = 12
  )

  # Overwrite defaults with any user-provided arguments in ...
  user_args <- list(...)
  final_args <- utils::modifyList(default_args, user_args)

  # Call $sample() using do.call
  fit <- do.call(mod$sample, final_args)


  # Model parameter summary
  all_pars<-c(
    "alphaCEonSur1Sur2" ,
    "b1CEonSur1Sur2",
    "b2CEonSur1Sur2",
    "SigSqCEonSur1Sur2",

    "alphaSur1onSur2",
    "bSur1onSur2",
    "SigSqSur1onSur2",

    "muSur2",
    "sigSqSur2",

    "R2CEonSur1Sur2", #R^2 for regression on CE~sur1 + sur2
    "RMSECEonSur1Sur2" #RMSE for regression on CE~sur1 + sur2
  )

  # Intersect with parameters present in this specific model
  #Exclude "alphaCEonSur1Sur2" if running fixed intercept model.
  present_pars <- intersect(all_pars, fit$metadata()$stan_variables)



  my_summary<-fit$summary(variables =  present_pars,
                          posterior::default_summary_measures()[1:4],
                          ~posterior::quantile2(.x, probs = c(0.025, 0.975))
                          )


  # 8. LOO and WAIC calculations

  # "log_lik" variable will be in the output.
  # You can now extract it directly.
  log_lik_matrix<- fit$draws("log_lik", format = "array")
  r_eff <- loo::relative_eff(
    x = exp(log_lik_matrix)
  )

  # efficient approximate leave-one-out (LOO) cross-validation for Bayesian models using Pareto smoothed importance sampling (PSIS)
  loo_result <- loo::loo(log_lik_matrix, r_eff = r_eff)

  #Widely applicable information criterion (WAIC)
  #compute WAIC from the pointwise log-likelihood
  waic_result<-loo::waic(log_lik_matrix)

  return(list(
    diagnostics=fit$diagnostic_summary(),
    posterior_samples=fit$draws(format = "df"),
    summary=my_summary,
    loo=loo_result,
    waic=waic_result,
    fit=fit
  ))


}
