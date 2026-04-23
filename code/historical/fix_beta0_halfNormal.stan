// fix beta0=0 with half normal

data {
  int<lower=1> N;                        // number of studies
  array[N] vector[3] obs_mean;           // observed treatment effects (CE, chronic, acute)
  array[N] matrix[3, 3] obs_var;         // study-specific observed covariance matrices
}
parameters{

  //Population mean
  real muSur2;              // mean for acute slope

  //meta-regression for modeling chronic slope ~ acute slope
  real alphaSur1onSur2;     // intercept
  real bSur1onSur2;         // slope

//meta-regression for modeling clinical endpoints ~ chronic slope + acute slope
  //real alphaCEonSur1Sur2;   // intercept (= beta0). Because it is treated as 0, we should just remove this.
  real b1CEonSur1Sur2;      // slope for chronic slope 
  real b2CEonSur1Sur2;      // slope for acute slope 
  
  // Standard deviation components
  real<lower=0> sigSur2;            //Standard deviation for acute slope
  real<lower=0> SigSur1onSur2;      //square root of conditional variance for modeling chronic slope
  real<lower=0> SigCEonSur1Sur2;    //square root of conditional variance for modeling clinical endpoint
  

  // standard normal distribution to construct psi later
  matrix[3, N] z;
}

transformed parameters {
  
  // Construct 'psi' here using the non-centered parameterization 
  array[N] vector[3] psi;  // Latent values per study (now an array of vectors)


  // Construct psi using the reparameterization formula: psi = mu + sqrt(variance) * z
  // We loop through each study (column of z)
  for (i in 1:N) {

    psi[i][3]=muSur2+sigSur2*z[3,i];
    psi[i][2]=alphaSur1onSur2 + bSur1onSur2*psi[i][3]+SigSur1onSur2*z[2,i];
    psi[i][1]= b1CEonSur1Sur2*  psi[i][2]   + b2CEonSur1Sur2* psi[i][3]+SigCEonSur1Sur2*z[1,i];
  }              
}




model {
  // Priors
  muSur2 ~ normal(0, 100);
  sigSur2~ normal(0, 5);
  
  alphaSur1onSur2 ~ normal(0, 100);
  bSur1onSur2 ~ normal(0, 100);
  SigSur1onSur2~ normal(0, 5);
  
  b1CEonSur1Sur2 ~ normal(0, 100);
  b2CEonSur1Sur2 ~ normal(0, 100);
  SigCEonSur1Sur2 ~ normal(0, 5);
  

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
    // The log-likelihood of the observed data. Used for leave-one-out cross-validation or other model comparison metrics.
    log_lik[i] = multi_normal_lpdf(obs_mean[i] | psi[i], obs_var[i]);
  }
  
  // Derive variances from the standard deviation parameters
  real sigSqSur2 = square(sigSur2);
  real SigSqSur1onSur2 = square(SigSur1onSur2);
  real SigSqCEonSur1Sur2 = square(SigCEonSur1Sur2);
}
