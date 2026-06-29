#This function is basically a wrapper + a bit more code
#That produces posterior statistics when we observe a new RCT
#We want posterior mean, SD, and percentiles.
#I will also produce samples, not just statistics.

#Different from v2:
#We allow to put 0 for intercept.


#Generate mean and covariance using prior data. We will use this to construct weights to update prior \xi.
MCMC_mean_cov_ChrAcu<-function(MCMC_dat,diffuse_se=100, diffuse=TRUE, beta0=TRUE){
  #MCMC_dat = data table
  #diffuse = logical. if we want to change the MCMC samples to have diffuse mean & Variance structure 
  #diffuse_se = SE of the variance for surrogates.
  #beta0 = logical. If TRUE, then the distribution of the intercept term of the meta-regression relating the treatment effect of CE to other surrogates will be centered to 
  
  library(data.table)
  library(mvtnorm)
  

  #Make sure we are using data.table format
  this.MCMC_dat<-data.table(MCMC_dat)
  
  #Save number of MCMC samples
  B<-nrow(this.MCMC_dat)
  
  
  #If we are specifying diffuse prior, set the corresponding MCMC parameters to what we want to set. 
  if(diffuse){
    this.MCMC_dat[,bSur1onSur2:=0]
    this.MCMC_dat[,SigSqSur1onSur2:=diffuse_se^2]
    this.MCMC_dat[,sigSqSur2:=diffuse_se^2]
    this.MCMC_dat[,muSur2:=0]
    this.MCMC_dat[,alphaSur1onSur2:=0]
  }
  
  #If we want to shift the mean of the intercept to 0
  if(beta0){
    this.MCMC_dat[,alphaCEonSur1Sur2:=scale(alphaCEonSur1Sur2, center=TRUE, scale=FALSE)]
  }
  
  #Variance 
  #define variance structure 
  this.MCMC_dat[,Sigma33:=sigSqSur2]
  this.MCMC_dat[,Sigma22:=bSur1onSur2^2 * sigSqSur2 + SigSqSur1onSur2]
  this.MCMC_dat[,Sigma23:=bSur1onSur2*sigSqSur2]
  this.MCMC_dat[,Sigma13:=b1CEonSur1Sur2*(bSur1onSur2*sigSqSur2) + b2CEonSur1Sur2*sigSqSur2]
  this.MCMC_dat[,Sigma12:=b1CEonSur1Sur2*(bSur1onSur2^2*sigSqSur2 + SigSqSur1onSur2)+ b2CEonSur1Sur2*bSur1onSur2*sigSqSur2]
  this.MCMC_dat[,Sigma11:=b1CEonSur1Sur2^2 *(bSur1onSur2^2*sigSqSur2 + SigSqSur1onSur2)+b2CEonSur1Sur2^2 * sigSqSur2 + SigSqCEonSur1Sur2 + 2*b1CEonSur1Sur2*b2CEonSur1Sur2*(bSur1onSur2*sigSqSur2)]
  
  
  
  #Array => list is fastest without using parallel
  
  
  
  #An array of covariances
  these_cols <- c("Sigma11","Sigma12","Sigma22","Sigma13","Sigma23","Sigma33")
  tmp <- as.matrix(this.MCMC_dat[, ..these_cols])   # B x 6
  covars_array <- array(0, dim = c(3, 3, B))
  # Fill upper triangle incl diag
  ut <- upper.tri(matrix(0, 3, 3), diag = TRUE)
  for (i in seq_len(B)) {
    covars_array[,,i][ut] <- tmp[i, ]
    covars_array[,,i][!ut] <- t(covars_array[,,i])[!ut]
  }
  #Convert to a list
  covars <- lapply(seq_len(dim(covars_array)[3]), function(i) covars_array[,,i])
  
  #Do likewise for means
  this.MCMC_dat[,mu1:=alphaCEonSur1Sur2 + b1CEonSur1Sur2*(alphaSur1onSur2 + bSur1onSur2*muSur2) + b2CEonSur1Sur2*muSur2]
  this.MCMC_dat[,mu2:=alphaSur1onSur2 + bSur1onSur2*muSur2]
  this.MCMC_dat[,mu3:=muSur2]
  tmp<-as.matrix(this.MCMC_dat[,.(mu1, mu2, mu3)])
  means_array <- array(0, dim = c(3, 1, B))
  for (i in seq_len(B)) {
    means_array[,,i] <- tmp[i, ]
  }
  means <- lapply(seq_len(dim(means_array)[3]), function(i) means_array[,,i,drop=FALSE])
  
 
  return(list(mean=means, covar=covars, dt=this.MCMC_dat[,.SD,.SDcols=c("mu1","mu2","mu3",these_cols)]))
  
  
}





theta0_given_hat_psi0<-function(mcmc_dat, sample_dat,diffuse=TRUE, diffuse_se=100, beta0=TRUE){
  #sample_dat = list of variable value.
  #mcmc_dat = data table
  #diffuse = logical. if we want to change the MCMC samples to have diffuse mean & Variance structure 
  #diffuse_se = SE of the variance for surrogates.
  #beta0 = logical. If TRUE, then the distribution of the intercept term of the meta-regression relating the treatment effect of CE to other surrogates will be centered to 
  
  library(parallel)
  
  this.MCMC_dat<-data.table(mcmc_dat)
  
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
  prior_mean_cov<-MCMC_mean_cov_ChrAcu(
    MCMC_dat = this.MCMC_dat,
    diffuse_se = diffuse_se,
    diffuse = diffuse,
    beta0 = beta0
  )

  
  #Calculate importance weight 
  w<-
    lapply(1:B,function(i){
      dmvnorm(x = hat_psi0, mean=prior_mean_cov$mean[[i]], sigma = hat_Sigma_y0 + prior_mean_cov$covar[[i]])
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
    rmvnorm(1, mean = post_mean_noW[[i]], sigma = post_var_noW[[i]])
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