#We will produce posterior samples from diffuse prior on surrogate
#Using surrogate MBE
#inv gamma
#Also, shift the distribution of beta0 to center at 0.
#Should set up 66 array
#And good number of cores per each



rm(list=ls())
library(data.table)
library(ggplot2)
library(dplyr)
library(parallel)
library(mvtnorm)
library(gt)
library(gtsummary)
library(gridExtra)
library(grid)
library(openxlsx)
library(tictoc)
ncores<-strtoi(Sys.getenv("SLURM_NTASKS")) #Pick up -ntasks or --n from the environment



# Setup job array in CHPC -------------------------------------------------

# Get command line arguments, expecting the task ID
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No task ID provided. Usage: Rscript prepare_data.R <TASK_ID>", call. = FALSE)
}

# Convert the argument to an integer row index
index <- as.integer(args[1])



# Data --------------------------------------------------------------------

#Read in all rds names
filenames<-list.files("../../MCMC_output_80000_repar/", pattern = ".rds")
#Exclude ALL_STUDIES
filenames<-filenames[!grepl("ALL_STUDIES", filenames)]
#Exclude the rds files with loo and waic in it
filenames<-filenames[!grepl("loo", filenames)]
filenames<-filenames[!grepl("waic", filenames)]


#choose the study
rds_file_name<-filenames[index]
study_name<-sub("*.rds$", "", rds_file_name) #Remove everything but the study name
study_name<-stringr::str_replace(study_name,"__","/")

#Read in the MCMC data
this.MCMC_dat<-readRDS(paste0("../../MCMC_output_80000_repar/",rds_file_name))
this.MCMC_dat<-data.table(this.MCMC_dat)

#Read in the study data that was left out
myData<-readRDS("../evt3data.rds") #to append the excluded data as the new RCT
sample_dat<-list(
  
  #Clinical endpoint
  ClnEst = myData[allrxname==study_name,loghr1_est],
  ClnSE = myData[allrxname==study_name,loghr1_se],
  
  #Chornic slope
  Sur1Est = myData[allrxname==study_name,beta31_beta32_est],
  Sur1SE = myData[allrxname==study_name,beta31_beta32_se],
  
  #acute slope
  Sur2Est = myData[allrxname==study_name,zbeta11_beta12_est],
  Sur2SE = myData[allrxname==study_name,zbeta11_beta12_se],
  
  R1Clin = myData[allrxname==study_name,beta3_ce_corr],
  R2Clin = myData[allrxname==study_name,beta1_ce_corr],
  R1R2 = myData[allrxname==study_name,dgfr_chr_corr]
)



# Run analysis ------------------------------------------------------------


#Math is derived in overleaf. Can be found in:
#CKD_two_surrogates/surrogate_MBE_math_Derivation.tex



  
diffuse_se<-100
this.MCMC_dat[,alphaCEonSur1Sur2:=scale(alphaCEonSur1Sur2,center=T,scale=F)]
  
#Calculate the values that are not relevant to the MCMC samples
  
#Generate the observed estimate and variance matrix.
hat_psi0<-matrix(c(sample_dat$Sur1Est,sample_dat$Sur2Est), nrow=2)
hat_Sigma_y0<-matrix(
  c(sample_dat$Sur1SE^2, sample_dat$Sur1SE*sample_dat$Sur2SE*sample_dat$R1R2,
    sample_dat$Sur1SE*sample_dat$Sur2SE*sample_dat$R1R2, sample_dat$Sur2SE^2),
  nrow=2,ncol = 2, byrow=T)
inv_hat_Sigma_y0<-solve(hat_Sigma_y0)


#Generate posterior distribution of surrogates conditional on new data \hat{\gamma}_0 
#Prior is bivariate standard normal distribution
#So the posterior is normal-normal conjugate update
post_gamma_var<-solve(diag(2)/diffuse_se^2 + inv_hat_Sigma_y0)
post_gamma_mean<-post_gamma_var %*% inv_hat_Sigma_y0 %*% hat_psi0
  
  
# first, let's sample the posterior gamma's from here 
# and get the posterior samples of \theta using those. 
post_gamm_samples<-rmvnorm(nrow(this.MCMC_dat), mean = post_gamma_mean, sigma = post_gamma_var)
post_ce_samples<-mclapply(1:nrow(this.MCMC_dat), function(i){
  this.MCMC_dat[
    i, 
    rnorm(1,
          mean=alphaCEonSur1Sur2+b1CEonSur1Sur2*post_gamm_samples[i,1]+b2CEonSur1Sur2*post_gamm_samples[i,2],
          sd=sqrt(SigSqCEonSur1Sur2))]
}, mc.cores = ncores)
post_ce_samples<-do.call(c,post_ce_samples)

#Estimate the marginal mean of the posterior \theta_{01}
#For each MCMC sample,
cond_means<-
  mclapply(1:nrow(this.MCMC_dat), function(i){
    this.MCMC_dat[i,alphaCEonSur1Sur2+b1CEonSur1Sur2*post_gamma_mean[1,1]+b2CEonSur1Sur2*post_gamma_mean[2,1]]
  }, mc.cores = ncores)
cond_means<-do.call(c,cond_means)
post_mean<-mean(cond_means) #this is the mean!


#Estimate the marginal variance of the posterior \theta_{01}
#Use Law of total variance

#Variance of conditional mean:
var_cond_mean<-var(cond_means) #Sample variance

#Expected value of of conditional variance:
mean_cond_var<-
  mclapply(1:nrow(this.MCMC_dat), function(i){
    this.MCMC_dat[i,b1CEonSur1Sur2^2*post_gamma_var[1,1]+b2CEonSur1Sur2^2*post_gamma_var[2,2]+2*b1CEonSur1Sur2*b2CEonSur1Sur2*post_gamma_var[1,2] + SigSqCEonSur1Sur2]
  }, mc.cores = ncores)

mean_cond_var<-Reduce("+", mean_cond_var) #Sum up the conditional variances
mean_cond_var<-mean_cond_var/nrow(this.MCMC_dat) #Average over all samples

#Sum them
post_var<-var_cond_mean + mean_cond_var

  # #Return mean and variance
  # output<-c(post_mean, post_var)
  # names(output)<-c("est","var")
  # return(output)
#   
#   return(
#     list(
#       "mean"=post_mean,
#       "var"=post_var,    
#       "post_samples"= post_ce_samples
# 
#       )
#   )
# })

# saveRDS(surr_MBE_out, file ="./surrogate_MBE_diffuse/inv_gamma/mean_var_samples.rds")
  
surr_MBE_out<-
  list(
    "mean"=post_mean,
    "var"=post_var,    
    "post_samples"= post_ce_samples
    
  )


#To save, replace / with __
save_name<-stringr::str_replace_all(study_name,"/","__")
saveRDS(surr_MBE_out, file = paste0( save_name, ".rds"))



