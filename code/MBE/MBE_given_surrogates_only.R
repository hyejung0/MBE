#Example generating distribution of $\theta_0 \mid \hat{\gamma}_0$
#under model D. 

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




