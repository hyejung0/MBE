#2025 Nov 05 (Wed)
#R code to run MCMC on all 66 EMA studies
#Shared in Tom's MBE statistics team
#surrogates: chronic slope & acute slope
#clinical endpoint: KFRT or GFR $ \leq 15 ml/min/1.73m^2$


rm(list=ls())
library(data.table)
library(tictoc)
# install.packages("remotes")
remotes::install_github("stan-dev/cmdstanr")


# Data --------------------------------------------------------------------

myData<-readRDS("evt3data.rds")




# Step 1: Stack observed means per study
obs_mean<-lapply(1:nrow(myData),function(j){
  unlist(myData[j,.(loghr1_est,            ## observed Clinical endpoint
                    beta31_beta32_est,     ## observed chronic slope
                    zbeta11_beta12_est     ## observed acute slope
  )])
})


# Step 2: Construct observed covariance matrix per study
obs_var <- lapply(1:nrow(myData), function(j) {
  # vector of SEs
  sample_dat<-myData[j,.(loghr1_se, #SE for Clinical endpoint
                         beta31_beta32_se, #SE for chronic slope
                         zbeta11_beta12_se, #SE for acute slope
                         beta3_ce_corr, #correlation of chronic slope & clinical endpoint
                         beta1_ce_corr, #correlation of acute slope & clinical endpoint
                         dgfr_chr_corr  #correlation of acute & chronic slope
  )]
  
  setnames(
    sample_dat,
    old=c("loghr1_se","beta31_beta32_se","zbeta11_beta12_se","beta3_ce_corr","beta1_ce_corr","dgfr_chr_corr"),
    new=c("ClnSE","Sur1SE","Sur2SE","R1Clin","R2Clin","R1R2")
  )
  
  
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

# define stan file --------------------------------------------------------

#Fixed prior for inverse gamma for variance structures. 
#You must change them manually if you want to modify them


stan_code_string<-"
data {
  int<lower=1> N;                        // number of studies
  array[N] vector[3] obs_mean;           // observed treatment effects (CE, chronic, acute)
  array[N] matrix[3, 3] obs_var;         // study-specific observed covariance matrices
}
parameters{

  //Population means
  real muSur2;
  real alphaSur1onSur2;
  real bSur1onSur2;
  real alphaCEonSur1Sur2;
  real b1CEonSur1Sur2;
  real b2CEonSur1Sur2;
  
  // variance components
  real<lower=0> sigSqSur2;
  real<lower=0> SigSqSur1onSur2;
  real<lower=0> SigSqCEonSur1Sur2;

  // *** CHANGE 1: Replace 'psi' with standard normal variates 'z' ***
  // Instead of sampling the correlated 'psi' directly, we sample 'z'
  // from a standard normal distribution. 'z' has 3 rows and N columns.
  matrix[3, N] z;
}

transformed parameters {
  
  // *** CHANGE 2: Construct 'psi' here using the non-centered parameterization ***
  array[N] vector[3] psi;  // Latent values per study (now an array of vectors)


  // Construct psi using the reparameterization formula: psi = mu + sqrt(variance) * z
  // We loop through each study (column of z)
  for (i in 1:N) {

    psi[i][3]=muSur2+sqrt(sigSqSur2)*z[3,i];
    psi[i][2]=alphaSur1onSur2 + bSur1onSur2*psi[i][3]+sqrt(SigSqSur1onSur2)*z[2,i];
    psi[i][1]=alphaCEonSur1Sur2 + b1CEonSur1Sur2*  psi[i][2]   + b2CEonSur1Sur2* psi[i][3]+sqrt(SigSqCEonSur1Sur2)*z[1,i];
  }              
}




model {
  // Priors
  muSur2 ~ normal(0, 100);
  sigSqSur2~ inv_gamma(0.261, 0.005);
  
  alphaSur1onSur2 ~ normal(0, 100);
  bSur1onSur2 ~ normal(0, 100);
  SigSqSur1onSur2~ inv_gamma(0.261, 0.005);
  
  alphaCEonSur1Sur2 ~ normal(0, 100);
  b1CEonSur1Sur2 ~ normal(0, 100);
  b2CEonSur1Sur2 ~ normal(0, 100);
  SigSqCEonSur1Sur2 ~ inv_gamma(0.261, 0.000408);
  

  to_vector(z) ~ std_normal();

  // The Stage 1 model (psi ~ multi_normal(mu, Sigma)) is now implicitly
  // defined by the construction of psi in the transformed parameters block.

  // Stage 2 model (Likelihood) - this remains the same
  for (i in 1:N) {
    obs_mean[i] ~ multi_normal(psi[i], obs_var[i]);
  }
  
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    // The log-likelihood of the observed data
    log_lik[i] = multi_normal_lpdf(obs_mean[i] | psi[i], obs_var[i]);
  }
}


"

# Write the Stan code string to a temporary file
# This returns a path to the temporary .stan file
temp_stan_file <- write_stan_file(code = stan_code_string)

# Compile the model from the temporary file
mod <- cmdstan_model(stan_file = temp_stan_file)




# Run --------------------------------------------------------------------

ncores<-3 #Number of cores to run chains parallely
options(mc.cores = ncores) #necessary for setting parallel 
num.iter=100000
num.warmup=num.iter/2
num.chains<-ncores
num.thin=1


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
  output_dir = "./run_MCMC_output/",  #This directory must exist first.
  output_basename="ALL_STUDIES" #the STAN files are saved with "ALL_STUDIES" as a base name
)
toc()


#Save posterior samples
draws_df <- fit$draws(format = "df")
saveRDS(draws_df,"./run_MCMC_output/ALL_STUDIES.rds") 



# Check convergence -------------------------------------------------------


#Re-build the STAN output from the saved csv files
output_files <- list.files("./run_MCMC_output", pattern = ".csv", full.names = T)
fit<-as_cmdstan_fit(output_files)


#Select the variables of interest. Define the variable names
all_pars<-c(
  
  #coefficients and variance for regression of CE ~ chronic slope + acute slope
  "alphaCEonSur1Sur2" ,
  "b1CEonSur1Sur2",
  "b2CEonSur1Sur2",  
  "SigSqCEonSur1Sur2",
  
  #coefficients and variance for regression of chronic slope ~ acute slope
  "alphaSur1onSur2", 
  "bSur1onSur2",       
  "SigSqSur1onSur2",   
  
  #marginal mean and variance for acute slope
  "muSur2",            
  "sigSqSur2"
  
  # ,"log_lik[1]" ,.. If you want to save log-likelihood, then save it.
)

#Print out the MCMC summary
my_summary<-fit$summary(variables = all_pars)
my_summary
# print(my_summary, n=10)


#Or we can add custom function to the summary:
my_summary<-fit$summary(variables =  all_pars,
                        posterior::default_summary_measures()[1:4],
                        quantiles = ~ quantile(., probs = c(0.025,0.05,0.1,0.9,0.95, 0.975)),
                        posterior::default_convergence_measures()
)
my_summary
