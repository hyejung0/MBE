#This model has highly coupled parameters and latent variables. Gibbs struggles when:
# Posteriors are strongly correlated
# Parameters are near non-identifiability
# You're already doing a great job approximating Stan's joint model, but Gibbs simply cannot explore high-dimensional coupled spaces efficiently unless all conditionals are tightly informative.

#Excellent — you're moving toward a powerful idea: block Metropolis-within-Gibbs is often used when:
#parameters are strongly correlated in the posterior, and
#full conditional distributions are not of standard form or are inefficient to sample separately.


#My regression coefficients are having problem. 
#Let's sample them together.
library(bayesplot)
library(data.table)
library(invgamma)
library(coda)
library(mvtnorm)

gibbs_sampler_v2 <- function(data, inits,burn.out=0.2,n_iter = 10000, parallel=FALSE,ncores=length(inits)) {
  #data = data frame containing the data
  #inits = list of paramter values to initalize. If more than one list if provided, then the number of inits list is number of chains. 
  #burn.out = proportion between 0 and 1, the proportion of n_iter to throw out as warm up.
  #parallel = logical. If true, then we will run parallel on chains
  #ncores= number of cores.
  
  m <- nrow(data)
  
  
  
  
  #For each observation, build a list of means and variances
  data_mean<-
    lapply(1:nrow(data),function(i){
      
      data[
        i,
        matrix(
          c(
            ClnEst,
            Sur1Est,
            Sur2Est
          ),
          nrow=3,ncol=1,byrow=T)
      ]
    })
  
  
  data_var<-
    lapply(1:nrow(data),function(i){
      
      data[
        i,
        matrix(
          c(ClnSE^2, R1Clin*ClnSE*Sur1SE,  R2Clin*ClnSE*Sur2SE, 
            R1Clin*ClnSE*Sur1SE, Sur1SE^2, Sur1SE*Sur2SE*R1R2,
            R2Clin*ClnSE*Sur2SE, Sur1SE*Sur2SE*R1R2,Sur2SE^2),
          nrow=3,ncol=3,byrow=T)
      ]
    })
  
  #calculate inverse
  inv_data_var<-lapply(data_var,solve)
  
  
  par_num<-length(inits[[1]]) + 2 #number of parameters, plus R2 and OptTotal, the generated parameters
  nchains<-length(inits)  #number of chains 
  
  #If we have more than one chain and parallel was specified as TRUE, run parallel
  if(parallel==TRUE & nchains>1){
    
    library(foreach)
    library(doParallel)
    
    
    # 1. Decide how many workers you want:
    n.cores <- min(ncores,nchains)
    
    # 2. Create a cluster.  Explicitly specify PSOCK for cross-platform compatibility.
    cl <- makeCluster(n.cores, type = "PSOCK")
    
    # 3. Register it so foreach knows to use it
    registerDoParallel(cl)
    
    
    # 4. Run your parallel loop!
    parallel_run <- foreach(this_chain = 1:nchains,
                            # how to stitch outputs together
                            .multicombine = FALSE,      # for modest # of cores you can leave FALSE
                            .packages = c("bayesplot", "data.table", "invgamma", "coda", "MASS")   # list any packages your workers need
    ) %dopar% {
      
      #for each chain...:
      
      
      # Storage
      draws <- list(
        muSur2 = numeric(n_iter),
        sigSqSur2 = numeric(n_iter),
        alphaSur1onSur2 = numeric(n_iter),
        bSur1onSur2 = numeric(n_iter),
        SigSqSur1onSur2 = numeric(n_iter),
        alphaCEonSur1Sur2 = numeric(n_iter),
        b1CEonSur1Sur2 = numeric(n_iter),
        b2CEonSur1Sur2 = numeric(n_iter),
        SigSqCEonSur1Sur2 = numeric(n_iter),
        R2 = numeric(n_iter),
        OptTotal = numeric(n_iter)
        
      )
      
      #Prior parmeter set up
      #For diffuse normal distribution
      prior_mean <- 0
      prior_var <- 100^2
      #For inverse gamma
      ig_shape <- 0.261
      ig_scale_logHR <- 0.000408
      ig_scale_GFR <- 0.005
      
      
      
      # Initialize
      list2env(inits[[this_chain]], envir = environment()) #this allows all elements in the list "sample_dat" to be available by its name.
      
      
      
      #Fun gibbs sampling
      for (iter in 1:n_iter) {
        
        
        
        #generate latent variables
        #This is just normal-normal conjugacy
        
        latent_mean<-c(
          #clinical endpoint
          alphaCEonSur1Sur2 + b1CEonSur1Sur2*(alphaSur1onSur2 + bSur1onSur2*muSur2) + b2CEonSur1Sur2*muSur2,
          
          # chronic slope
          alphaSur1onSur2 + bSur1onSur2*muSur2,
          
          #acute slope
          muSur2
          
        )
        
        #Latent covariance
        Sigma_latent <- matrix(0, 3, 3)
        Sigma_latent[3,3]<-sigSqSur2
        Sigma_latent[2,2] <-  bSur1onSur2^2 * sigSqSur2 + SigSqSur1onSur2
        Sigma_latent[2,3]<-Sigma_latent[3,2] <- bSur1onSur2*sigSqSur2
        Sigma_latent[1,3]<-Sigma_latent[3,1] <- b1CEonSur1Sur2*(bSur1onSur2*sigSqSur2) + b2CEonSur1Sur2*sigSqSur2
        Sigma_latent[1,2]<-Sigma_latent[2,1] <-b1CEonSur1Sur2*(bSur1onSur2^2*sigSqSur2 + SigSqSur1onSur2)+ b2CEonSur1Sur2*bSur1onSur2*sigSqSur2
        Sigma_latent[1,1]<-b1CEonSur1Sur2^2 *(bSur1onSur2^2*sigSqSur2 + SigSqSur1onSur2)+b2CEonSur1Sur2^2 * sigSqSur2 + SigSqCEonSur1Sur2 + 2*b1CEonSur1Sur2*b2CEonSur1Sur2*(bSur1onSur2*sigSqSur2)
        inv_Sigma_latent<- solve(Sigma_latent)
        
        #Update the latent variable using normal-normal conjugacy
        post_covar<-lapply(inv_data_var, function(xx)solve(xx+inv_Sigma_latent))
        
        #Posterior mean
        post_mu<-
          lapply(1:length(post_covar),function(i){
            post_covar[[i]]%*%inv_data_var[[i]]%*%data_mean[[i]] +
              post_covar[[i]]%*%inv_Sigma_latent%*%latent_mean
          })
        
        # Sample from posterior
        psi_i <- 
          sapply(1:length(post_mu),function(i){
            MASS::mvrnorm(1, post_mu[[i]], post_covar[[i]])
          })
        # Extract observed values
        ClnEst <- psi_i[1,]
        Sur1Est <- psi_i[2,]
        Sur2Est <- psi_i[3,]
        
        #Use these latent variables to update the hyperparameters
        
        
        ### 1. muSur2 | rest
        sigma_post <- 1 / ( (m / sigSqSur2) + (1 / prior_var))
        mu_post <- sigma_post * (sum(Sur2Est) / sigSqSur2)
        muSur2 <- rnorm(1, mu_post, sqrt(sigma_post))
        
        ### 2. sigSqSur2 | rest
        shape_post <- ig_shape + m / 2
        resids<-Sur2Est - muSur2
        scale_post <- ig_scale_GFR + (sum(resids^2) / 2)
        sigSqSur2 <- rinvgamma(1, shape_post, rate = scale_post)
        
        
        
        ### 3, 4. alphaSur1onSur2 & bSur1onSur2| rest
        X<-cbind(1,Sur2Est)
        covar_mat<-solve(diag(ncol(X))/prior_var + t(X)%*%X/SigSqSur1onSur2)
        mean_mat<-covar_mat %*% (t(X) %*%matrix(Sur1Est,ncol=1))/SigSqSur1onSur2
        sample_together<-mvtnorm::rmvnorm(n=1,mean=mean_mat,sigma = covar_mat)
        alphaSur1onSur2<-sample_together[1,1]
        bSur1onSur2<-sample_together[1,2]
        
        
        
        # ### 3. alphaSur1onSur2 | rest
        # sigma_post <- 1 / (m / SigSqSur1onSur2 + 1 / prior_var)
        # c_i<-Sur1Est - bSur1onSur2*Sur2Est
        # mu_post <- sigma_post * (sum(c_i) / SigSqSur1onSur2)
        # alphaSur1onSur2 <- rnorm(1, mu_post, sqrt(sigma_post))
        # 
        # ### 4. bSur1onSur2 | rest 
        # sigma_post <- 1 / ( (sum(Sur2Est^2) / SigSqSur1onSur2) + (1 / prior_var))
        # b_c_i<-Sur2Est*(Sur1Est - alphaSur1onSur2)
        # mu_post <- sigma_post * (sum(b_c_i) / SigSqSur1onSur2)
        # bSur1onSur2 <- rnorm(1, mu_post, sqrt(sigma_post))
        
        
        ### 5. SigSqSur1onSur2 | rest
        shape_post <- ig_shape +( m / 2)
        resid<-Sur1Est - alphaSur1onSur2 - bSur1onSur2*Sur2Est
        scale_post <- ig_scale_GFR + (sum(resid^2) / 2)
        SigSqSur1onSur2 <- rinvgamma(1, shape_post, rate = scale_post)

        
        ### 6.7.8. alphaCEonSur1Sur2 & b1CEonSur1Sur2 & b2CEonSur1Sur2 | rest
        X<-cbind(1,Sur1Est,Sur2Est)
        covar_mat<-solve(diag(ncol(X))/prior_var + t(X)%*%X/SigSqCEonSur1Sur2)
        mean_mat<-covar_mat %*% (t(X) %*%matrix(ClnEst,ncol=1))/SigSqCEonSur1Sur2
        sample_together<-mvtnorm::rmvnorm(n=1,mean=mean_mat,sigma = covar_mat)
        alphaCEonSur1Sur2<-sample_together[1,1]
        b1CEonSur1Sur2<-sample_together[1,2]
        b2CEonSur1Sur2<-sample_together[1,3]
        
        # ### 6. alphaCEonSur1Sur2 | rest
        # sigma_post <- 1 / ((m / SigSqCEonSur1Sur2) + (1 / prior_var))
        # b_c_i<-ClnEst - b1CEonSur1Sur2*Sur1Est - b2CEonSur1Sur2*Sur2Est
        # mu_post <- sigma_post * (sum(b_c_i) / SigSqCEonSur1Sur2)
        # alphaCEonSur1Sur2 <- rnorm(1, mu_post, sqrt(sigma_post))
        # 
        # ### 7. b1CEonSur1Sur2 | rest
        # sigma_post <- 1 / (( sum(Sur1Est^2) / SigSqCEonSur1Sur2) + (1 / prior_var))
        # b_c_i<-Sur1Est*(ClnEst - alphaCEonSur1Sur2 - b2CEonSur1Sur2*Sur2Est)
        # mu_post <- sigma_post * (sum(b_c_i) / SigSqCEonSur1Sur2)
        # b1CEonSur1Sur2 <- rnorm(1, mu_post, sqrt(sigma_post))
        # 
        # ### 8. b2CEonSur1Sur2 | rest
        # sigma_post <- 1 / (( sum(Sur2Est^2) / SigSqCEonSur1Sur2) + (1 / prior_var))
        # b_c_i<-Sur2Est*(ClnEst - alphaCEonSur1Sur2 - b1CEonSur1Sur2*Sur1Est)
        # mu_post <- sigma_post * (sum(b_c_i) / SigSqCEonSur1Sur2)
        # b2CEonSur1Sur2 <- rnorm(1, mu_post, sqrt(sigma_post))
        
        
        ### 9. SigSqCEonSur1Sur2 | rest
        shape_post <- ig_shape + m / 2
        resid <- ClnEst - (alphaCEonSur1Sur2 + b1CEonSur1Sur2 * Sur1Est +  b2CEonSur1Sur2 * Sur2Est)
        scale_post <- ig_scale_logHR + ( sum(resid^2) / 2)
        SigSqCEonSur1Sur2 <- rinvgamma(1, shape_post, rate = scale_post)
        
        
        
        
        # resids_mat[iter,]<-Sur2Est - muSur2
        
        ### Generated paramters
        
        Var_clin<- 
          b1CEonSur1Sur2^2 *(bSur1onSur2^2*sigSqSur2 + SigSqSur1onSur2)+ 
          b2CEonSur1Sur2^2 * sigSqSur2 + SigSqCEonSur1Sur2 + 
          2*b1CEonSur1Sur2*b2CEonSur1Sur2*(bSur1onSur2*sigSqSur2)
        
        R2<-1 - SigSqCEonSur1Sur2/Var_clin
        
        OptTotal<- ((b1CEonSur1Sur2/b2CEonSur1Sur2) * 3 + 3)/12
        
        
        
        
        # Save draws
        draws$muSur2[iter] <- muSur2
        draws$sigSqSur2[iter] <- sigSqSur2
        draws$alphaSur1onSur2[iter] <- alphaSur1onSur2
        draws$bSur1onSur2[iter] <- bSur1onSur2
        draws$SigSqSur1onSur2[iter] <- SigSqSur1onSur2
        draws$alphaCEonSur1Sur2[iter] <- alphaCEonSur1Sur2
        draws$b1CEonSur1Sur2[iter] <- b1CEonSur1Sur2
        draws$b2CEonSur1Sur2[iter] <- b2CEonSur1Sur2
        draws$SigSqCEonSur1Sur2[iter] <- SigSqCEonSur1Sur2
        draws$R2[iter]<-R2
        draws$OptTotal[iter]<-OptTotal
      }  
      
      draws_table<-as.data.table(draws)
      
      
      #calculate statistics
      
      #burn.out?
      if(burn.out!=0){
        draws_table<-draws_table[-c(1:nrow(draws_table)*burn.out),]
      }
      
      #return the data
      list(
        "bind.draws"=draws_table,
        "effective.size"=effectiveSize(draws_table)
      )
      
    }
    
    # 5. Clean up
    stopCluster(cl)
    
    
    
    #For each parallel run, find the same object:
    
    #posterior samples
    bind.draws<-lapply(parallel_run,"[[","bind.draws")
    bind.draws<-rbindlist(bind.draws)
    
    #effective sample size
    effective.size<-lapply(parallel_run,"[[","effective.size")
    effective.size<-do.call(rbind,effective.size)
    
    
  }else{
    #If we don't want to do parallel, don't do parallel.
    bind.draws<-NULL
    effective.size<-matrix(0,nrow=nchains,ncol=par_num)
    
    
    for(this_chain in 1:nchains){#for each chain:
      
      
      # Storage
      draws <- list(
        muSur2 = numeric(n_iter),
        sigSqSur2 = numeric(n_iter),
        alphaSur1onSur2 = numeric(n_iter),
        bSur1onSur2 = numeric(n_iter),
        SigSqSur1onSur2 = numeric(n_iter),
        alphaCEonSur1Sur2 = numeric(n_iter),
        b1CEonSur1Sur2 = numeric(n_iter),
        b2CEonSur1Sur2 = numeric(n_iter),
        SigSqCEonSur1Sur2 = numeric(n_iter),
        R2 = numeric(n_iter),
        OptTotal = numeric(n_iter)
        
      )
      
      
      #Prior parmeter set up
      #For diffuse normal distribution
      prior_mean <- 0
      prior_var <- 100^2
      #For inverse gamma
      ig_shape <- 0.261
      ig_scale_small <- 0.000408
      ig_scale_medium <- 0.005
      
      
      
      # Initialize
      list2env(inits[[this_chain]], envir = environment()) #this allows all elements in the list "sample_dat" to be available by its name.
      
      
      # resids_mat<-matrix(NA, nrow=n_iter, ncol=m)
      for (iter in 1:n_iter) {
        
        
        
        #generate latent variables
        #This is just normal-normal conjugacy
        
        latent_mean<-c(
          #clinical endpoint
          alphaCEonSur1Sur2 + b1CEonSur1Sur2*(alphaSur1onSur2 + bSur1onSur2*muSur2) + b2CEonSur1Sur2*muSur2,
          
          # chronic slope
          alphaSur1onSur2 + bSur1onSur2*muSur2,
          
          #acute slope
          muSur2
          
        )
        
        #Latent covariance
        Sigma_latent <- matrix(0, 3, 3)
        Sigma_latent[3,3]<-sigSqSur2
        Sigma_latent[2,2] <-  bSur1onSur2^2 * sigSqSur2 + SigSqSur1onSur2
        Sigma_latent[2,3]<-Sigma_latent[3,2] <- bSur1onSur2*sigSqSur2
        Sigma_latent[1,3]<-Sigma_latent[3,1] <- b1CEonSur1Sur2*(bSur1onSur2*sigSqSur2) + b2CEonSur1Sur2*sigSqSur2
        Sigma_latent[1,2]<-Sigma_latent[2,1] <-b1CEonSur1Sur2*(bSur1onSur2^2*sigSqSur2 + SigSqSur1onSur2)+ b2CEonSur1Sur2*bSur1onSur2*sigSqSur2
        Sigma_latent[1,1]<-b1CEonSur1Sur2^2 *(bSur1onSur2^2*sigSqSur2 + SigSqSur1onSur2)+b2CEonSur1Sur2^2 * sigSqSur2 + SigSqCEonSur1Sur2 + 2*b1CEonSur1Sur2*b2CEonSur1Sur2*(bSur1onSur2*sigSqSur2)
        inv_Sigma_latent<- solve(Sigma_latent)
        
        #Update the latent variable using normal-normal conjugacy
        post_covar<-lapply(inv_data_var, function(xx)solve(xx+inv_Sigma_latent))
        
        #Posterior mean
        post_mu<-
          lapply(1:length(post_covar),function(i){
            post_covar[[i]]%*%inv_data_var[[i]]%*%data_mean[[i]] +
              post_covar[[i]]%*%inv_Sigma_latent%*%latent_mean
          })
        
        # Sample from posterior
        psi_i <- 
          sapply(1:length(post_mu),function(i){
            MASS::mvrnorm(1, post_mu[[i]], post_covar[[i]])
          })
        # Extract observed values
        ClnEst <- psi_i[1,]
        Sur1Est <- psi_i[2,]
        Sur2Est <- psi_i[3,]
        
        #Use these latent variables to update the hyperparameters
        
        
        ### 1. muSur2 | rest
        sigma_post <- 1 / ( (m / sigSqSur2) + (1 / prior_var))
        mu_post <- sigma_post * (sum(Sur2Est) / sigSqSur2)
        muSur2 <- rnorm(1, mu_post, sqrt(sigma_post))
        
        ### 2. sigSqSur2 | rest
        shape_post <- ig_shape + m / 2
        resids<-Sur2Est - muSur2
        scale_post <- ig_scale_small + (sum(resids^2) / 2)
        sigSqSur2 <- rinvgamma(1, shape_post, rate = scale_post)
        
        
        ### 3, 4. alphaSur1onSur2 & bSur1onSur2| rest
        X<-cbind(1,Sur2Est)
        covar_mat<-solve(diag(ncol(X))/prior_var + t(X)%*%X/SigSqSur1onSur2)
        mean_mat<-covar_mat %*% (t(X) %*%matrix(Sur1Est,ncol=1))/SigSqSur1onSur2
        sample_together<-mvtnorm::rmvnorm(n=1,mean=mean_mat,sigma = covar_mat)
        alphaSur1onSur2<-sample_together[1,1]
        bSur1onSur2<-sample_together[1,2]

        
        
        # ### 3. alphaSur1onSur2 | rest
        # sigma_post <- 1 / (m / SigSqSur1onSur2 + 1 / prior_var)
        # c_i<-Sur1Est - bSur1onSur2*Sur2Est
        # mu_post <- sigma_post * (sum(c_i) / SigSqSur1onSur2)
        # alphaSur1onSur2 <- rnorm(1, mu_post, sqrt(sigma_post))
        # 
        # ### 4. bSur1onSur2 | rest 
        # sigma_post <- 1 / ( (sum(Sur2Est^2) / SigSqSur1onSur2) + (1 / prior_var))
        # b_c_i<-Sur2Est*(Sur1Est - alphaSur1onSur2)
        # mu_post <- sigma_post * (sum(b_c_i) / SigSqSur1onSur2)
        # bSur1onSur2 <- rnorm(1, mu_post, sqrt(sigma_post))
        
        
        ### 5. SigSqSur1onSur2 | rest
        shape_post <- ig_shape +( m / 2)
        resid<-Sur1Est - alphaSur1onSur2 - bSur1onSur2*Sur2Est
        scale_post <- ig_scale_medium + (sum(resid^2) / 2)
        SigSqSur1onSur2 <- rinvgamma(1, shape_post, rate = scale_post)
        
        
        
        ### 6.7.8. alphaCEonSur1Sur2 & b1CEonSur1Sur2 & b2CEonSur1Sur2 | rest
        X<-cbind(1,Sur1Est,Sur2Est)
        covar_mat<-solve(diag(ncol(X))/prior_var + t(X)%*%X/SigSqCEonSur1Sur2)
        mean_mat<-covar_mat %*% (t(X) %*%matrix(ClnEst,ncol=1))/SigSqCEonSur1Sur2
        sample_together<-mvtnorm::rmvnorm(n=1,mean=mean_mat,sigma = covar_mat)
        alphaCEonSur1Sur2<-sample_together[1,1]
        b1CEonSur1Sur2<-sample_together[1,2]
        b2CEonSur1Sur2<-sample_together[1,3]
        
        
        # ### 6. alphaCEonSur1Sur2 | rest
        # sigma_post <- 1 / ((m / SigSqCEonSur1Sur2) + (1 / prior_var))
        # b_c_i<-ClnEst - b1CEonSur1Sur2*Sur1Est - b2CEonSur1Sur2*Sur2Est
        # mu_post <- sigma_post * (sum(b_c_i) / SigSqCEonSur1Sur2)
        # alphaCEonSur1Sur2 <- rnorm(1, mu_post, sqrt(sigma_post))
        # 
        # ### 7. b1CEonSur1Sur2 | rest
        # sigma_post <- 1 / (( sum(Sur1Est^2) / SigSqCEonSur1Sur2) + (1 / prior_var))
        # b_c_i<-Sur1Est*(ClnEst - alphaCEonSur1Sur2 - b2CEonSur1Sur2*Sur2Est)
        # mu_post <- sigma_post * (sum(b_c_i) / SigSqCEonSur1Sur2)
        # b1CEonSur1Sur2 <- rnorm(1, mu_post, sqrt(sigma_post))
        # 
        # ### 8. b2CEonSur1Sur2 | rest
        # sigma_post <- 1 / (( sum(Sur2Est^2) / SigSqCEonSur1Sur2) + (1 / prior_var))
        # b_c_i<-Sur2Est*(ClnEst - alphaCEonSur1Sur2 - b1CEonSur1Sur2*Sur1Est)
        # mu_post <- sigma_post * (sum(b_c_i) / SigSqCEonSur1Sur2)
        # b2CEonSur1Sur2 <- rnorm(1, mu_post, sqrt(sigma_post))
        
        
        ### 9. SigSqCEonSur1Sur2 | rest
        shape_post <- ig_shape + m / 2
        resid <- ClnEst - (alphaCEonSur1Sur2 + b1CEonSur1Sur2 * Sur1Est +  b2CEonSur1Sur2 * Sur2Est)
        scale_post <- ig_scale_small + ( sum(resid^2) / 2)
        SigSqCEonSur1Sur2 <- rinvgamma(1, shape_post, rate = scale_post)
        
        
        
        
        # resids_mat[iter,]<-Sur2Est - muSur2
        
        ### Generated paramters
        
        Var_clin<- 
          b1CEonSur1Sur2^2 *(bSur1onSur2^2*sigSqSur2 + SigSqSur1onSur2)+ 
          b2CEonSur1Sur2^2 * sigSqSur2 + SigSqCEonSur1Sur2 + 
          2*b1CEonSur1Sur2*b2CEonSur1Sur2*(bSur1onSur2*sigSqSur2)
        
        R2<-1 - SigSqCEonSur1Sur2/Var_clin
        
        OptTotal<- ((b1CEonSur1Sur2/b2CEonSur1Sur2) * 3 + 3)/12
        
        
        
        
        # Save draws
        draws$muSur2[iter] <- muSur2
        draws$sigSqSur2[iter] <- sigSqSur2
        draws$alphaSur1onSur2[iter] <- alphaSur1onSur2
        draws$bSur1onSur2[iter] <- bSur1onSur2
        draws$SigSqSur1onSur2[iter] <- SigSqSur1onSur2
        draws$alphaCEonSur1Sur2[iter] <- alphaCEonSur1Sur2
        draws$b1CEonSur1Sur2[iter] <- b1CEonSur1Sur2
        draws$b2CEonSur1Sur2[iter] <- b2CEonSur1Sur2
        draws$SigSqCEonSur1Sur2[iter] <- SigSqCEonSur1Sur2
        draws$R2[iter]<-R2
        draws$OptTotal[iter]<-OptTotal
      }  
      
      
      draws_table<-as.data.table(draws)
      
      
      #burn.out?
      if(burn.out!=0){
        draws_table<-draws_table[-c(1:nrow(draws_table)*burn.out),]
      }
      
      bind.draws<-rbind(bind.draws,draws_table)
      effective.size[this_chain,]<-effectiveSize(draws_table)
      
      
    }
    
    colnames(effective.size)<-colnames(draws_table)
    
  }
  
  
  
  
  
  #Calculate all calculation..
  effective.sample.size<-apply(effective.size,2,sum)
  quantile_par<-apply(bind.draws,2,function(x)quantile(x,c(.025,.25,.5,.75,.975)))
  quantile_par<-t(quantile_par)
  
  
  
  
  #calculate the estimated marginal posterior variance of the estimand.
  #Get the number of iteration and number of chains.
  
  #number of iteration:
  n<-(n_iter - n_iter*burn.out)/2
  
  #number of chains
  m<-nchains*2
  
  #Bind it to the sample:
  bind.draws[,chains:=rep(1:m, each=n)]
  
  #for each chain, calculate sample mean and variance
  all.sample.mean<-bind.draws[,lapply(.SD,mean),by=chains]
  all.sample.mean[,chains:=NULL]
  all.sample.var<-bind.draws[,lapply(.SD,var),by=chains]
  all.sample.var[,chains:=NULL]
  
  
  Rhat<-sapply(1:par_num,function(i){
    sample.mean<-all.sample.mean[[i]]
    sample.var<-all.sample.var[[i]]
    
    #Between variance
    B.par<-n/(m-1)*sum((sample.mean-mean(sample.mean))^2)
    
    #Within variance
    W.par<-mean(sample.var)
    
    #Total variance
    var.hat<-((n-1)/n)*W.par + B.par/n
    
    #Return RHat
    sqrt(var.hat[1]/W.par)
  })
  names(Rhat)<-colnames(all.sample.mean)
  
  
  
  list("effective.size"=effective.sample.size, 
       "Rhat"=Rhat,
       "quantile.par"=quantile_par,
       "sample.mean"=all.sample.mean,
       "sample.var"=all.sample.var,
       "draws"=bind.draws)
  
  
}
