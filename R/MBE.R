#' Estimate an MBE Posterior with Two Surrogate Endpoints
#'
#' @description Uses new-trial estimates and historical posterior draws to
#' estimate the posterior distribution of an MBE by importance sampling. The
#' model is fixed to one clinical endpoint and exactly two surrogate endpoints.
#'
#' @param mcmc_dat A data frame or matrix containing posterior draws from the
#'   two-surrogate historical model.
#' @param sample_dat A named list containing exactly nine scalar values:
#'   `ClnEst`, `ClnSE`, `Sur1Est`, `Sur1SE`, `Sur2Est`, `Sur2SE`, `R1Clin`,
#'   `R2Clin`, and `R12`.
#' @param diffuse_se A positive numeric value specifying the standard deviation
#'   of the diffuse prior for the two surrogate effects.
#' @param diffuse A logical value indicating whether to use a diffuse prior. Defaults to TRUE.
#' @param intercept0 A logical value indicating whether to center the intercept term of the meta-regression to zero. Defaults to TRUE.
#'
#' @return A list containing posterior means, posterior covariance, weighted
#'   quantiles for the clinical and two surrogate effects, posterior draws, the
#'   corresponding importance weights, and importance-sampling diagnostics.
#'   The diagnostic vector reports the effective sample size, effective sample
#'   size relative to the number of historical draws, and maximum normalized
#'   weight. Endpoint components are always ordered as clinical endpoint,
#'   surrogate 1, and surrogate 2.
#'
#' @export
#'
#' @examples
#' # Use the example data provided in the package.
#' # Load the historical posterior samples where the model was fit
#' # with random intercept and inverse gamma prior for variance parameters
#' # on the 65/66 simulated trials in `trial_sim_dat` dataset.
#' # It's the first trial that was left out as if it was a new trial.
#' # We demonstrate estimating the posterior distribution of MBE on this first
#' # row that was left out.
#' data("historical_posterior")
#' data("trial_sim_dat")
#' one_sim_dat<-list(
#'
#' # treatment effect on CE
#' ClnEst=trial_sim_dat$CE_est[1],
#' ClnSE=trial_sim_dat$CE_se[1],
#'
#' # treatment effect on chronic slope
#' Sur1Est=trial_sim_dat$Sur1_est[1],
#' Sur1SE=trial_sim_dat$Sur1_se[1],
#'
#' # treatment effect on acute slope
#' Sur2Est=trial_sim_dat$Sur2_est[1],
#' Sur2SE=trial_sim_dat$Sur2_se[1],
#'
#' # correlation between CE and chronic slope
#' R1Clin=trial_sim_dat$Cor_CE_Sur1[1],
#'
#' # correlation between CE and acute slope
#' R2Clin=trial_sim_dat$Cor_CE_Sur2[1],
#'
#' # correlation between chronic slope and acute slope
#' R12=trial_sim_dat$Cor_Sur1_Sur2[1]
#' )
#'
#' MBE_distribution<-MBE(
#' mcmc_dat = historical_posterior,
#' sample_dat=one_sim_dat,
#' diffuse_se = 100,
#' diffuse = TRUE,
#' intercept0 = TRUE
#' )
#'
#' head(MBE_distribution$post_psi0)
MBE <- function(mcmc_dat, sample_dat, diffuse = TRUE, diffuse_se = 100, intercept0 = TRUE) {
  .validate_mcmc_dat(mcmc_dat)
  .validate_sample_dat(sample_dat)
  .validate_scalar_logical(diffuse, "diffuse")
  .validate_scalar_logical(intercept0, "intercept0")
  .validate_positive_scalar(diffuse_se, "diffuse_se")

  this.MCMC_dat <- data.table::as.data.table(data.table::copy(mcmc_dat))

  #save the total number of MCMC samples
  B <- nrow(this.MCMC_dat)


  #Construct observed mean and variance vector and matrix:
  hat_psi0 <- c(sample_dat$ClnEst, sample_dat$Sur1Est, sample_dat$Sur2Est)
  hat_Sigma_y0 <- matrix(
    c(sample_dat$ClnSE^2, sample_dat$R1Clin*sample_dat$ClnSE*sample_dat$Sur1SE,  sample_dat$R2Clin*sample_dat$ClnSE*sample_dat$Sur2SE,
      sample_dat$R1Clin*sample_dat$ClnSE*sample_dat$Sur1SE, sample_dat$Sur1SE^2, sample_dat$Sur1SE*sample_dat$Sur2SE*sample_dat$R12,
      sample_dat$R2Clin*sample_dat$ClnSE*sample_dat$Sur2SE, sample_dat$Sur1SE*sample_dat$Sur2SE*sample_dat$R12,sample_dat$Sur2SE^2),
    nrow = 3, ncol = 3, byrow = TRUE
  )
  .validate_positive_definite(hat_Sigma_y0, "The sampling covariance matrix derived from `sample_dat`")
  inv_hat_Sigma_y0 <- chol2inv(chol(hat_Sigma_y0))

  #Let's first re-parameterize the MCMC parameters into marginal mean and variance.
  prior_mean_cov<-vector_to_matrix(
    MCMC_dat = this.MCMC_dat,
    diffuse_se = diffuse_se,
    diffuse = diffuse,
    intercept0 = intercept0
  )


  #Calculate importance weight
  log_w <- vapply(seq_len(B), function(i) {
    mvtnorm::dmvnorm(
      x = hat_psi0,
      mean = prior_mean_cov$mean[[i]],
      sigma = hat_Sigma_y0 + prior_mean_cov$covar[[i]],
      log = TRUE
    )
  }, numeric(1))
  if (anyNA(log_w) || !any(is.finite(log_w))) {
    stop(
      "The importance weights could not be normalized; check the posterior draws and new-trial inputs.",
      call. = FALSE
    )
  }
  w <- exp(log_w - max(log_w))
  norm_w <- w / sum(w)
  importance_ess <- 1 / sum(norm_w^2)
  importance_diagnostics <- c(
    effective_sample_size = importance_ess,
    relative_effective_sample_size = importance_ess / B,
    maximum_normalized_weight = max(norm_w)
  )




  #draw a sample of \psi_0 from conditional posterior distribution, \pi(\psi_0 \mid \hat{\psi}_0, \; \xi(1)^{(b)}) <- incorporating new observation
  #This can be simply calculated using normal-normal conjugacy.


  #construct variance:

  #First, take inverse of the prior mean
  inv_prior_cov <- lapply(prior_mean_cov$covar, function(xx) {
    chol2inv(chol(xx))
  })
  post_var_noW <- lapply(inv_prior_cov, function(xx) {
    chol2inv(chol(xx + inv_hat_Sigma_y0))
  })


  #construct mean
  post_mean_noW<-lapply(1:B,function(i){
    post_var_noW[[i]] %*%
      (inv_hat_Sigma_y0 %*% hat_psi0 + inv_prior_cov[[i]] %*% prior_mean_cov$mean[[i]])
  })


  #draw a sample using the mean and variance.
  almost_post_psi0<-lapply(1:B,function(i){
    mvtnorm::rmvnorm(1, mean = post_mean_noW[[i]], sigma = post_var_noW[[i]])
  })



  #generate posterior samples by re-drawing sample with normalzed weight as the probability.
  #select sample index
  sample_idx<-sample(1:B, size = B, replace = TRUE, prob = norm_w)
  post_psi0<-almost_post_psi0[sample_idx]
  post_psi0 <- do.call(rbind, post_psi0)
  colnames(post_psi0) <- c("clinical", "surrogate1", "surrogate2")


  #calculate posterior statistics using weighted samples technique


  #Posterior mean
  post_mean<-lapply(1:B,function(i){
    post_mean_noW[[i]] * norm_w[i] #multiply with weight
  })
  post_mean <- Reduce("+", post_mean) #weighted sum


  #Posterior variance is to be calculated in two parts:
  #1. expected conditional variance
  #2. variance of conditional means

  #1. expected conditional variance
  exp_cond_var<-lapply(1:B,function(i){
    post_var_noW[[i]] * norm_w[i] #multiply with weight
  })
  exp_cond_var <- Reduce("+", exp_cond_var) #weighted sum

  #2. variance of conditional means
  var_cond_mean<-lapply(1:B,function(i){
    norm_w[i] * ( (post_mean_noW[[i]] - post_mean) %*% t(post_mean_noW[[i]] - post_mean) )
  })
  var_cond_mean <- Reduce("+", var_cond_mean) #weighted sum

  #3. Total variance
  post_var <- exp_cond_var + var_cond_mean
  endpoint_names <- c("clinical", "surrogate1", "surrogate2")
  rownames(post_mean) <- endpoint_names
  dimnames(post_var) <- list(endpoint_names, endpoint_names)


  #Calculate quantiles
  #1. Sort the samples in the order smallest to largest
  psi0<-do.call(rbind,almost_post_psi0)
  colnames(psi0) <- c("clinical", "surrogate1", "surrogate2")
  psi0 <- data.table::as.data.table(psi0)
  psi0[,norm_w:=norm_w]
  psi0[,w:=w]
  psi01 <- psi0[, .(clinical, norm_w)][order(clinical)]
  psi02 <- psi0[, .(surrogate1, norm_w)][order(surrogate1)]
  psi03 <- psi0[, .(surrogate2, norm_w)][order(surrogate2)]

  #2. empirical CDF is cumulative sum of the weights
  psi01[,empirical_cdf:=cumsum(norm_w)]
  psi02[,empirical_cdf:=cumsum(norm_w)]
  psi03[,empirical_cdf:=cumsum(norm_w)]

  #3.Find the first sample that is greater than or equivalent to quantile of interest.
  these_quant<-
    c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)
  quant01<-sapply(these_quant,function(this.quantile){
    psi01[empirical_cdf >= this.quantile, clinical][1]
  })
  quant02<-sapply(these_quant,function(this.quantile){
    psi02[empirical_cdf >= this.quantile, surrogate1][1]
  })
  quant03<-sapply(these_quant,function(this.quantile){
    psi03[empirical_cdf >= this.quantile, surrogate2][1]
  })
  names(quant01)<-
    names(quant02)<-
    names(quant03)<-
    paste0("quantile_",these_quant)




  list(
    post_mean = post_mean,
    post_var = post_var,
    post_quantiles_clinical = quant01,
    post_quantiles_surrogate1 = quant02,
    post_quantiles_surrogate2 = quant03,
    post_psi0 = post_psi0,
    weight_data = psi0,
    importance_diagnostics = importance_diagnostics
  )
}
