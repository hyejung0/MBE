#Calculate MBE
#Should have samples saved from historical MCMC run

rm(list=ls())
library(data.table)


# Data --------------------------------------------------------------------

#Read in the MCMC data
this.MCMC_dat<-readRDS("./run_MCMC_output/ALL_STUDIES.rds")
this.MCMC_dat<-data.table(this.MCMC_dat)

#Read in the data from the new RCT that the MBE will be applied to
myData<-readRDS("evt3data.rds") # This is a data.table
study_name<-"CANVAS" #suppose this is the study we want to use as a new study


#generate list of all variables needed for MBE calculation
sample_dat<-list(
  
  #Clinical endpoint
  ClnEst = myData[allrxname==study_name,loghr1_est], #CE estimate
  ClnSE = myData[allrxname==study_name,loghr1_se], #CE standard error
  
  #Chronic slope
  Sur1Est = myData[allrxname==study_name,beta31_beta32_est], #chronic slope estimate
  Sur1SE = myData[allrxname==study_name,beta31_beta32_se], #chronic slope standard error
  
  #acute slope
  Sur2Est = myData[allrxname==study_name,zbeta11_beta12_est], #acute slope estimate
  Sur2SE = myData[allrxname==study_name,zbeta11_beta12_se],  #acute slope standard error
  
  #Correlations
  R1Clin = myData[allrxname==study_name,beta3_ce_corr], #CE and chronic slope
  R2Clin = myData[allrxname==study_name,beta1_ce_corr], #CE and acute slope
  R1R2 = myData[allrxname==study_name,dgfr_chr_corr] #chronic & acute slope
)






# Calculate MBE -----------------------------------------------------------

#read in the function
source("post_theta0_samples_v4.R")


MBE_full_post<-post_theta0_samples_v4(
  mcmc_dat = this.MCMC_dat, #the MCMC data set from historical data fitting
  sample_dat = sample_dat, #The new trial data
  diffuse=TRUE, #Logical. If we want to put diffuse prior on the surrogate endpoints, set TRUE. FALSE otherwise.
  diffuse_se = 100, #If diffuse=TRUE, we should specify the SE of the surrogate endpoints for the diffuse normal prior.
  intercept0=FALSE #Logical. If we want to set the mean of the intercept of the meta-regression of treatment effect on CE on the surrogate endpoint, set TRUE.
  ) 

#This function returns:
# post_mean  = posterior mean (mathematically calculated)
# post_var = posterior variance (mathematically calculated)
# post_quantiles_psi01 = percentiles of posterior distribution on CE
# post_quantiles_psi02 = percentiles of posterior distribution on chronic slope
# post_quantiles_psi03 = percentiles of posterior distribution on acute slope
# post_psi0 = posterior distribution sample (weighted)
# weight_data = data.table containing unweighted posterior draws and the associated weights for those draws. Please use the norm_w (normalized weight). 


MBE_full_post$post_var


MBE_full_post$weight_data

MBE_full_post$post_mean
MBE_full_post$post_var

