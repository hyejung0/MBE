#' Estimate posterior distribution of Multi-Component Bayesian Endpoints (MBE)
#'
#' @description Takes a numeric vector and returns a data frame containing
#' the mean, standard deviation, and count of non-missing values.
#'
#' @param mcmc_dat A data.table object containing MCMC samples.
#' @param sample_dat A list containing the observed mean and variance vector.
#' @param diffuse_se A numeric value specifying the standard error for the diffuse prior.
#' @param diffuse A logical value indicating whether to use a diffuse prior. Defaults to TRUE.
#' @param intercept0 A logical value indicating whether to center the intercept term of the meta-regression to zero. Defaults to TRUE.
#' stripped before the computation proceeds. Defaults to `TRUE`.
#' @param ncores An integer. For parallel computing.
#'
#' @return A data frame with three columns: `mean`, `sd`, and `n`.
#'
#' @export
#'
#' @examples
#' # Use the example MCMC data provided in the package.
#' # Use a single time point in the provided simulated data.
#'
#' # We will use simulation case 2, month 30.
#' sim_dat<-data.table::data.table(sim_dat)
#' sim_dat<-sim_dat[case==2 & analysis.month==30, ]
#' one_sim_dat<-list(
#'
#' # treatment effect on CE
#' ClnEst=sim_dat$ClnEst,
#' ClnSE=sim_dat$ClnSE,
#'
#' # treatment effect on chronic slope
#' Sur1Est=sim_dat$Sur1Est,
#' Sur1SE=sim_dat$Sur1SE,
#'
#' # treatment effect on acute slope
#' Sur2Est=sim_dat$Sur2Est,
#' Sur2SE=sim_dat$Sur2SE,
#'
#' # correlation between CE and chronic slope
#' R1Clin=sim_dat$R1Clin,
#'
#' # correlation between CE and acute slope
#' R2Clin=sim_dat$R2Clin,
#'
#' # correlation between chronic slope and acute slope
#' R1R2=sim_dat$R12
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
#'head(MBE_distribution)



MBE<-function(mcmc_dat, sample_dat,diffuse=TRUE, diffuse_se=100, ncores=3, intercept0=TRUE){


  this.MCMC_dat<-data.table::data.table(mcmc_dat)

  #save the total number of MCMC samples
  B<-nrow(this.MCMC_dat)


  #Construct observed mean and variance vector and matrix:
  hat_psi0<-c(sample_dat$ClnEst,
              sample_dat$Sur1Est,
              sample_dat$Sur2Est)
  hat_Sigma_y0<-matrix(
    c(sample_dat$ClnSE^2, sample_dat$R1Clin*sample_dat$ClnSE*sample_dat$Sur1SE,  sample_dat$R2Clin*sample_dat$ClnSE*sample_dat$Sur2SE,
      sample_dat$R1Clin*sample_dat$ClnSE*sample_dat$Sur1SE, sample_dat$Sur1SE^2, sample_dat$Sur1SE*sample_dat$Sur2SE*sample_dat$R1R2,
      sample_dat$R2Clin*sample_dat$ClnSE*sample_dat$Sur2SE, sample_dat$Sur1SE*sample_dat$Sur2SE*sample_dat$R1R2,sample_dat$Sur2SE^2),
    nrow=3,ncol = 3, byrow=T)
  inv_hat_Sigma_y0<-solve(hat_Sigma_y0)

  #Let's first re-parameterize the MCMC parameters into marginal mean and variance.
  prior_mean_cov<-vector_to_matrix(
    MCMC_dat = this.MCMC_dat,
    diffuse_se = diffuse_se,
    diffuse = diffuse,
    intercept0 = intercept0
  )


  #Calculate importance weight
  w<-
    lapply(1:B,function(i){
      mvtnorm::dmvnorm(x = hat_psi0, mean=prior_mean_cov$mean[[i]], sigma = hat_Sigma_y0 + prior_mean_cov$covar[[i]])
    })
  w<-do.call(c,w)
  #normalize the weight
  norm_w<-w/sum(w)




  #draw a sample of \psi_0 from conditional posterior distribution, \pi(\psi_0 \mid \hat{\psi}_0, \; \xi(1)^{(b)}) <- incorporating new observation
  #This can be simply calculated using normal-normal conjugacy.


  #construct variance:

  #First, take inverse of the prior mean
  inv_prior_cov<-lapply(prior_mean_cov$covar,function(xx){
    solve(xx)
  })
  post_var_noW<- lapply(inv_prior_cov,function(xx){
    solve(xx + inv_hat_Sigma_y0)
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
  post_psi0<-do.call(rbind,post_psi0)
  colnames(post_psi0)<-c("psi01","psi02","psi03")


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
  post_var<-exp_cond_var+var_cond_mean


  #Calculate quantiles
  #1. Sort the samples in the order smallest to largest
  psi0<-do.call(rbind,almost_post_psi0)
  colnames(psi0)<-c("psi01","psi02","psi03")
  psi0<-data.table(psi0)
  psi0[,norm_w:=norm_w]
  psi0[,w:=w]
  psi01<-psi0[,.(psi01,norm_w)][order(psi01,decreasing = F)]
  psi02<-psi0[,.(psi02,norm_w)][order(psi02,decreasing = F)]
  psi03<-psi0[,.(psi03,norm_w)][order(psi03,decreasing = F)]

  #2. empirical CDF is cumulative sum of the weights
  psi01[,empirical_cdf:=cumsum(norm_w)]
  psi02[,empirical_cdf:=cumsum(norm_w)]
  psi03[,empirical_cdf:=cumsum(norm_w)]

  #3.Find the first sample that is greater than or equivalent to quantile of interest.
  these_quant<-
    c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)
  quant01<-sapply(these_quant,function(this.quantile){
    psi01[empirical_cdf>=this.quantile,][1,psi01]
  })
  quant02<-sapply(these_quant,function(this.quantile){
    psi02[empirical_cdf>=this.quantile,][1,psi02]
  })
  quant03<-sapply(these_quant,function(this.quantile){
    psi03[empirical_cdf>=this.quantile,][1,psi03]
  })
  names(quant01)<-
    names(quant02)<-
    names(quant03)<-
    paste0("quantile_",these_quant)




  return(
    list(
      post_mean=post_mean, #posterior mean (mathematically calculated)
      post_var=post_var, #posterior variance (mathematically calculated)
      post_quantiles_psi01 = quant01, #percentiles of posterior distribution on CE
      post_quantiles_psi02 = quant02, #percentiles of posterior distribution on chronic slope
      post_quantiles_psi03 = quant03, #percentiles of posterior distribution on acute slope
      post_psi0=post_psi0, #posterior distribution sample (weighted)
      weight_data=psi0 #data.table containing unweighted posterior draws and the associated weights for those draws. Please use the norm_w (normalized weight).

    )
  )
}
