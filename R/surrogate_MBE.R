#' Estimate the MBE posterior distribution CE with only the two surrogate endpoints
#'
#' @description
#' Uses the estimated treatment effects on the two surrogate endpoints of
#' a new trial and the historical posterior draws to estimate the posterior
#' distribution of MBE.
#'
#' @param mcmc_dat A data frame or matrix containing posterior draws from the
#'   two-surrogate historical model. Must contains columns of the variables named:
#'   `alphaCEonSur1Sur2`, `b1CEonSur1Sur2`, `b2CEonSur1Sur2`, and `SigSqCEonSur1Sur2`.
#'
#' @param sample_dat A named list containing exactly nine scalar values:
#'   `Sur1Est`, `Sur1SE`, `Sur2Est`, `Sur2SE`, and `R12`.
#' @param diffuse A logical value indicating whether to use a diffuse prior. Defaults to TRUE.
#' @param diffuse_se A positive numeric value specifying the standard deviation
#'   of the diffuse prior for the two surrogate effects. Only used when `diffuse = TRUE`.
#'   defaults to 100.
#' @param intercept0 A logical value indicating whether to center the intercept
#' term of the meta-regression to zero. Defaults to TRUE.
#'
#'
#' @details
#' Let  \eqn{\hat{\gamma}_0 = (\hat{\gamma}_{0, 1}, \hat{\gamma}_{0, 2})^T} be the vector of
#' treatment effects on the two surrogate endpoints of a new trial, and let
#' \eqn{\hat{\gamma}_0 \sim N(\gamma_0, \Sigma_0)} where \eqn{\Sigma_0} is
#' assumed to be known (estimated from the new trial data). We model prior on
#' \eqn{\gamma_0 = (\gamma_{0, 1}, \gamma_{0, 2})^T}, the vector of true
#' treatment effects on the two surrogate endpoints. The posterior distribution
#' of \eqn{\gamma_0} is normal by conjugacy. The prior distribution can either be
#' estimated from the historical posterior draws or be diffuse. For diffuse, a bivariate
#' normal distribution with mean vector 0 and covariance matrix
#' `diffuse_se`^2 \eqn{ \times I_2} is used. The bigger the `diffuse_se`,
#' the more diffuse the prior is. Samples of the posterior distribution of \eqn{\gamma_0}
#' are obtained and then used to estimate the posterior distribution of MBE on the clinical endpoint
#' using the relationship:
#' \deqn{\theta_0 = \alpha_\theta + \beta_{\gamma_1} \cdot \gamma_{0, 1} +
#' \beta_{\gamma_2} \cdot \gamma_{0, 2} + \epsilon,
#' \quad \epsilon \sim N(0, \lambda^2_\theta).}
#'
#' Historical posterior distribution of \eqn{\alpha_\theta, \beta_{\gamma_1}, \beta_{\gamma_2}, \lambda^2_\theta}
#' is used to estimate the posterior distribution of MBE on the clinical endpoint without updating its distribution.

#' @return A data table containing:
#'   - `post_Sigma22`: Posterior draws of the variance of the first surrogate.
#'   - `post_Sigma23`: Posterior draws of the covariance between the first and second surrogates.
#'   - `post_Sigma33`: Posterior draws of the variance of the second surrogate.
#'   - `post_mu2`: Posterior draws of the mean of the first surrogate.
#'   - `post_mu3`: Posterior draws of the mean of the second surrogate.
#'   - `post_CE`: Posterior draws of the treatment effect on the clinical endpoint.
#'   - `post_psi2`: Posterior draws of the treatment effect on the first surrogate.
#'   - `post_psi3`: Posterior draws of the treatment effect on the second surrogate.
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

#' # treatment effect on chronic slope
#' Sur1Est=trial_sim_dat$Sur1_est[1],
#' Sur1SE=trial_sim_dat$Sur1_se[1],
#'
#' # treatment effect on acute slope
#' Sur2Est=trial_sim_dat$Sur2_est[1],
#' Sur2SE=trial_sim_dat$Sur2_se[1],
#'
#' # correlation between chronic slope and acute slope
#' R12=trial_sim_dat$Cor_Sur1_Sur2[1]
#' )
#'
#'MBE_from_surrogate_only<-surrogate_MBE(
#' mcmc_dat = historical_posterior,
#' sample_dat=one_sim_dat,
#' diffuse = TRUE, #We will not use the historical posterior distribution of the two surrogate endpoints to estimate the posterior distribution of MBE. Instead, we will use a diffuse prior.
#' diffuse_se = 100,
#' intercept0 = FALSE)
#'
surrogate_MBE<-function(mcmc_dat, sample_dat, diffuse=TRUE, diffuse_se=100, intercept0=TRUE){


  # Check if sample_dat has the required elements
  required_elements <- c("Sur1Est", "Sur1SE", "Sur2Est", "Sur2SE", "R12")
  if (!all(required_elements %in% names(sample_dat))) {
    stop("sample_dat must contain the following elements: ", paste(required_elements, collapse = ", "))
  }


  #convert mcmc_dat to data.table. Leave the original mcmc_dat. Just create a new one.
  this.MCMC_dat <- copy(data.table::data.table(mcmc_dat))

  #Check if this.MCMC_dat contains all required columns.
  if(diffuse){
    required_columns <- c("alphaCEonSur1Sur2", "b1CEonSur1Sur2", "b2CEonSur1Sur2", "SigSqCEonSur1Sur2")
  } else {
    required_columns <- c("alphaCEonSur1Sur2", "b1CEonSur1Sur2", "b2CEonSur1Sur2", "SigSqCEonSur1Sur2",
                          "alphaSur1onSur2", "bSur1onSur2", "SigSqSur1onSur2",
                          "muSur2","sigSqSur2")
  }
  if (!all(required_columns %in% colnames(this.MCMC_dat))) {
    stop("this.MCMC_dat must contain the following columns: ", paste(required_columns, collapse = ", "))
  }

  #If intercept0 is TRUE, center the intercept term of the meta-regression to zero
  if(intercept0){
    this.MCMC_dat[,alphaCEonSur1Sur2:=scale(alphaCEonSur1Sur2,center=T,scale=F)]
  }

  #If we are specifying diffuse prior, set the corresponding MCMC parameters to what we want to set.
  if(diffuse){
    this.MCMC_dat[,bSur1onSur2:=0]
    this.MCMC_dat[,SigSqSur1onSur2:=diffuse_se^2]
    this.MCMC_dat[,sigSqSur2:=diffuse_se^2]
    this.MCMC_dat[,muSur2:=0]
    this.MCMC_dat[,alphaSur1onSur2:=0]
  }

  #Number of MCMC samples
  B<-nrow(this.MCMC_dat)


  ###  Prior  ###
  #define variance structure
  this.MCMC_dat[,Sigma33:=sigSqSur2]
  this.MCMC_dat[,Sigma22:=bSur1onSur2^2 * sigSqSur2 + SigSqSur1onSur2]
  this.MCMC_dat[,Sigma23:=bSur1onSur2*sigSqSur2]

  #generate inverse of the 2 by 2 matrix of the two surrogate endpoints
  this.MCMC_dat[,det_Sigma22_33:=Sigma22*Sigma33 - Sigma23^2]
  this.MCMC_dat[,inv_Sigma22:=Sigma33/det_Sigma22_33]
  this.MCMC_dat[,inv_Sigma23:=-Sigma23/det_Sigma22_33]
  this.MCMC_dat[,inv_Sigma33:=Sigma22/det_Sigma22_33]

  #Do likewise for means
  this.MCMC_dat[,mu2:=alphaSur1onSur2 + bSur1onSur2*muSur2]
  this.MCMC_dat[,mu3:=muSur2]




  ###   Data   ###
  #Generate colmns of the observed data required to calculate posterior distribution
  #to this.MCMC_dat.
  this.MCMC_dat[,hat_psi02:=sample_dat$Sur1Est]
  this.MCMC_dat[,hat_psi03:=sample_dat$Sur2Est]
  det_Sigma_y0<-sample_dat$Sur1SE^2*sample_dat$Sur2SE^2 - (sample_dat$Sur1SE*sample_dat$Sur2SE*sample_dat$R12)^2
  this.MCMC_dat[,hat_Sigma_y0_22:=sample_dat$Sur2SE^2/det_Sigma_y0]
  this.MCMC_dat[,hat_Sigma_y0_23:=-(sample_dat$Sur1SE*sample_dat$Sur2SE*sample_dat$R12)/det_Sigma_y0]
  this.MCMC_dat[,hat_Sigma_y0_33:=sample_dat$Sur1SE^2/det_Sigma_y0]




  ### Posterior ###
  post_var<-data.table(
    entry22=this.MCMC_dat[,inv_Sigma22 + hat_Sigma_y0_22],
    entry23=this.MCMC_dat[,inv_Sigma23 + hat_Sigma_y0_23],
    entry33=this.MCMC_dat[,inv_Sigma33 + hat_Sigma_y0_33]
  ) #I need to take inverse of them
  post_var[,det:=entry22*entry33 - entry23^2]
  post_var[,inv_entry22:=entry33/det]
  post_var[,inv_entry23:=-entry23/det]
  post_var[,inv_entry33:=entry22/det]
  #Append the posteiror vairance to this.MCMC_dat
  this.MCMC_dat[,post_Sigma22:=post_var$inv_entry22]
  this.MCMC_dat[,post_Sigma23:=post_var$inv_entry23]
  this.MCMC_dat[,post_Sigma33:=post_var$inv_entry33]

  #Posterior mean
  a<-this.MCMC_dat[,hat_Sigma_y0_22*hat_psi02 + hat_Sigma_y0_23*hat_psi03]
  b<-this.MCMC_dat[,inv_Sigma22*mu2 + inv_Sigma23*mu3]

  c<-this.MCMC_dat[,hat_Sigma_y0_23*hat_psi02 + hat_Sigma_y0_33*hat_psi03]
  d<-this.MCMC_dat[,inv_Sigma23*mu2 + inv_Sigma33*mu3]

  almost_post_mu2<-a+b
  almost_post_mu3<-c+d

  this.MCMC_dat$post_mu2<-
  this.MCMC_dat$post_Sigma22*almost_post_mu2 + this.MCMC_dat$post_Sigma23*almost_post_mu3

  this.MCMC_dat$post_mu3<-
    this.MCMC_dat$post_Sigma23*almost_post_mu2 + this.MCMC_dat$post_Sigma33*almost_post_mu3


  #Using the posterior mean and variance of the surrogate endpoints, draw a new sample for each row
  this.MCMC_dat[
    ,
    c("post_psi2", "post_psi3") := as.list(
      mvtnorm::rmvnorm(
        1,
        mean = c(post_mu2, post_mu3),
        sigma = matrix(
          c(post_Sigma22, post_Sigma23,
            post_Sigma23, post_Sigma33),
          nrow = 2
        )
      )[1, ]
    ),
    by = .I
  ]


  #Generate posterior distribution of CE conditional on the posterior samples of the two surrogate endpoints
  this.MCMC_dat[
    ,
    post_CE :=
      alphaCEonSur1Sur2 +
      b1CEonSur1Sur2 * post_psi2 +
      b2CEonSur1Sur2 * post_psi3 +
      rnorm(.N, mean = 0, sd = sqrt(SigSqCEonSur1Sur2))
  ]

  #Return the posterior data only
  return(
    this.MCMC_dat[
      ,
      .(post_Sigma22, post_Sigma23, post_Sigma33, post_mu2, post_mu3, post_CE, post_psi2, post_psi3)
    ]
  )
}
