#Example code to generate 
#the distribution of $\theta_0 \mid \hat{\psi}_0$ under different model specifications
#When we incorporate all observed endpoints in the new RCT.
#Uses the historical posterior distribution.




rm(list=ls())
library(data.table)


setwd("/uufs/chpc.utah.edu/common/home/greene-group4/MBE/code/MBE")
source("theta0_given_hat_psi0.R") #Function to estimate the distribution of $\theta_0 \mid \hat{\psi}_0$

# Data --------------------------------------------------------------------


#Read in the historical posterior distribution of parameters 
historical_posteriors<-readRDS(
  "/uufs/chpc.utah.edu/common/home/greene-group4/Hyejung/MCMC_EVt3_ChronAcut/MCMC_output_80000_repar/Kamper.rds"
)
historical_posteriors<-data.table(historical_posteriors)



#Read in the new data. 
#save its data as a list of vectors. 

#Read in the study data that was left out
myData<-readRDS("../evt3data.rds") #this is 1 by 9 data table, with the 9 columns being the 9 elements as shown in 
#Data Structure section in README file.

sample_dat<-as.list(myData)





# Model A: carry over full historical posterior ---------------------------

modelA<-theta0_given_hat_psi0(
  mcmc_dat = historical_posteriors, #historical posterior samples 
  sample_dat = sample_dat, #List of data from the new RCT
  diffuse=FALSE, #Indicate whether or not we want to put diffuse MBE prior on the surrogate endpoints
  # diffuse_se = 100, #If diffuse=TRUE, this is the standard error of the diffuse MBE prior on the surrogate endpoints. The larger this number is, the more diffuse the prior is. This corresponds to the fixed value "c" in the variance in the manuscript This argument is ignored if diffuse=FALSE.
  beta0=FALSE) #Indicate whether or not we want to shift the mean of the beta the intercept to be zero. If FALSE, the historical posterior beta0 samples will be used as MBE's prior. If TRUE, the historical posterior beta0 samples will be centered to have sample mean zero before it is used as MBE's prior.



# Model B: Diffuse surrogate endpoints -------------------------------------

modelB<-theta0_given_hat_psi0(
  mcmc_dat = historical_posteriors, 
  sample_dat = sample_dat,
  diffuse = TRUE,
  diffuse_se = 100, # In the manuscript, c=100.
  beta0=FALSE
)



# Model C: beta0 fixed at zero --------------------------------------------

#This is the model where 
#the historical posterior was fitted with beta0 fixed at zero.
#We also use diffuse prior on the surrogate endpoint like in model D. 
modelC<-theta0_given_hat_psi0(
  mcmc_dat = historical_posteriors, 
  sample_dat = sample_dat,
  diffuse = TRUE,
  diffuse_se = 100, 
  beta0=FALSE
)



# Model D: shift mean of beta0 --------------------------------------------

#This one takes Model D but shifts the sample mean of beta0 to be zero.

modelD<-theta0_given_hat_psi0(
  mcmc_dat = historical_posteriors, 
  sample_dat = sample_dat,
  diffuse = TRUE,
  diffuse_se = 100,
  beta0=TRUE
)
