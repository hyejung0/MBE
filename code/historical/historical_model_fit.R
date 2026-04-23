
#R code to run MCMC on historical data
#surrogates: chronic slope & acute slope
#clinical endpoint: KFRT or GFR $ \leq 15 ml/min/1.73m^2$


rm(list=ls())
library(data.table)
library(tictoc)
# install.packages("remotes")
#remotes::install_github("stan-dev/cmdstanr")


# Data --------------------------------------------------------------------


#Read in the historical RCT data.
#The data structure follows the description in the README file.
myData<-readRDS("evt3data.rds")




# Step 1: Stack observed treatment effect per study
obs_mean<-lapply(1:nrow(myData),function(j){
  unlist(myData[j,.(ClnEst,      ## observed Clinical endpoint
                    Sur1Est,     ## observed chronic slope
                    Sur2Est      ## observed acute slope
  )])
})


# Step 2: Construct observed covariance matrix per study
obs_var <- lapply(1:nrow(myData), function(j) {
  # vector of SEs
  sample_dat<-myData[j,.(ClnSE, #SE for Clinical endpoint
                         Sur1SE, #SE for chronic slope
                         Sur2SE, #SE for acute slope
                         R1Clin, #correlation of chronic slope & clinical endpoint
                         R2Clin, #correlation of acute slope & clinical endpoint
                         R1R2  #correlation of acute & chronic slope
  )]

  
  sample_dat[
    ,
    matrix(
      c(ClnSE^2, R1Clin*ClnSE*Sur1SE,  R2Clin*ClnSE*Sur2SE, 
        R1Clin*ClnSE*Sur1SE, Sur1SE^2, Sur1SE*Sur2SE*R1R2,
        R2Clin*ClnSE*Sur2SE, Sur1SE*Sur2SE*R1R2,Sur2SE^2),
      nrow=3,ncol = 3, byrow=T)
  ]
})


#collect them together as a more comprehensive list
data_sub = list(
  N=nrow(myData), # set number of studies. should be 66.
  obs_mean = obs_mean,
  obs_var = obs_var
  
)  


# Statistical Modeling using STAN -----------------------------------------

#Source the stan file that you want to use. That is, run one of the four options below:

# Option 1 : fit the historical model with inverse gamma prior for variance parameters and $\beta_0$ fixed to be zero
mod <- cmdstan_model(stan_file = "fix_beta0_invGamma.stan")


#Option 2 : fit the historical model with inverse gamma prior for variance parameters and model $\beta_0$ as a random parameter
mod <- cmdstan_model(stan_file = "random_beta0_invGamma.stan")


# Option 3 : fit the historical model with half-normal prior for standard deviation parameters and $\beta_0$ fixed to be zero
mod <- cmdstan_model(stan_file = "fix_beta0_halfNormal.stan")

# Option 4 : fit the historical model with half-normal prior for standard deviation parameters and model $\beta_0$ as a random parameter
mod <- cmdstan_model(stan_file = "random_beta0_halfNormal.stan")




# Run --------------------------------------------------------------------

ncores<-3 #Number of cores to run chains parallely
options(mc.cores = ncores) #necessary for setting parallel 
num.iter=40000 #number of iterations for each chain. You can change this number as you like. We recommend to have at least 40,000 iterations for the historical model fitting.
num.warmup=num.iter/2 #number of warmup iterations. The default is half of the total number of iterations, but you can change this number as you like. We recommend to have at least 20,000 warmup iterations for the historical model fitting.
num.chains<-ncores
num.thin=1 #thinning number. The default is 1, which means no thinning. You can change this number as you like. We recommend to have no thinning for the historical model fitting.


set.seed(1) #Change seed number as you like
tic()
#Fit MCMC
fit <- mod$sample(
  data = data_sub,
  chains = num.chains,
  parallel_chains = num.chains,
  iter_warmup = num.warmup,
  iter_sampling = num.warmup,
  adapt_delta = 0.99,
  max_treedepth=25,
  output_dir = ".",  #current working directory is where the STAN output will be saved.
  output_basename="historical_fit" #the STAN files are saved with "historical_fit" as a base name
)
toc()


#Save posterior samples
draws_df <- fit$draws(format = "df")
saveRDS(draws_df,"./historical_fit.rds") 

