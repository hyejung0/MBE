data{

int nStudies;

vector[nStudies] chrn_eff;
vector[nStudies] chrn_se;
vector[nStudies] actslp_eff;
vector[nStudies] actslp_se;
vector[nStudies] ClnEst;
vector[nStudies] ClnSE;
// vector[nStudies] acr_eff;
// vector[nStudies] acr_se;

vector[nStudies] Rclnchron;
vector[nStudies] Rclnacute;
// vector[nStudies] Rclnacr;
vector[nStudies] Rchronacute;
// vector[nStudies] Rchronacr;
// vector[nStudies] Racuteacr;

}

parameters{

matrix[nStudies,3] withinStudyMeans;

real muSur2;
real<lower=0> sigSqSur2;

real alphaSur1onSur2;
real bSur1onSur2;
real<lower=0> SigSqSur1onSur2;

real alphaCEonSur1Sur2;
real b1CEonSur1Sur2;
real b2CEonSur1Sur2;
real<lower=0> SigSqCEonSur1Sur2;

}

transformed parameters {
  real<lower=0> sigSur2;
  real<lower=0> SigSur1onSur2;
  real<lower=0> SigCEonSur1Sur2;
  

  sigSur2=sqrt(sigSqSur2);
  SigSur1onSur2=sqrt(SigSqSur1onSur2);
  SigCEonSur1Sur2=sqrt(SigSqCEonSur1Sur2);
}

model{

for(jj in 1:nStudies){

vector[3] myY;
matrix[3,3] myVar;

myY[1]=ClnEst[jj];
myY[2]=chrn_eff[jj];
myY[3]=actslp_eff[jj];
// myY[4]=acr_eff[jj];

myVar[1,1]=ClnSE[jj]^2;
myVar[2,2]=chrn_se[jj]^2;
myVar[3,3]=actslp_se[jj]^2;
// myVar[4,4]=acr_se[jj]^2;

myVar[1,2]=Rclnchron[jj]*ClnSE[jj]*chrn_se[jj];
myVar[1,3]=Rclnacute[jj]*ClnSE[jj]*actslp_se[jj];
// myVar[1,4]=Rclnacr[jj]*ClnSE[jj]*acr_se[jj];
myVar[2,1]=Rclnchron[jj]*ClnSE[jj]*chrn_se[jj];
myVar[3,1]=Rclnacute[jj]*ClnSE[jj]*actslp_se[jj];
// myVar[4,1]=Rclnacr[jj]*ClnSE[jj]*acr_se[jj];

myVar[2,3]=Rchronacute[jj]*chrn_se[jj]*actslp_se[jj];
// myVar[2,4]=Rchronacr[jj]*chrn_se[jj]*acr_se[jj];
myVar[3,2]=Rchronacute[jj]*chrn_se[jj]*actslp_se[jj];
// myVar[4,2]=Rchronacr[jj]*chrn_se[jj]*acr_se[jj];

// myVar[3,4]=Racuteacr[jj]*actslp_se[jj]*acr_se[jj];
// myVar[4,3]=Racuteacr[jj]*actslp_se[jj]*acr_se[jj];

myY~multi_normal(withinStudyMeans[jj,],myVar);
                     
withinStudyMeans[jj,3]~normal(muSur2,sigSur2);

withinStudyMeans[jj,2]~normal(alphaSur1onSur2 + bSur1onSur2*withinStudyMeans[jj,3],SigSur1onSur2);

withinStudyMeans[jj,1]~normal(alphaCEonSur1Sur2 + b1CEonSur1Sur2*withinStudyMeans[jj,2] + b2CEonSur1Sur2*withinStudyMeans[jj,3],SigCEonSur1Sur2);

}



muSur2~normal(0,100);
sigSqSur2~inv_gamma(0.261,0.005);

alphaSur1onSur2~normal(0,100);
bSur1onSur2~normal(0,100);
SigSqSur1onSur2~inv_gamma(0.261,0.005);

alphaCEonSur1Sur2~normal(0,100);
b1CEonSur1Sur2~normal(0,100);
b2CEonSur1Sur2~normal(0,100);
SigSqCEonSur1Sur2~inv_gamma(0.261,0.000408);

}





// return R2 and OptTotal as well 
generated quantities {
  real<lower=0> Var_clin;
  real R2;
  real OptTotal;

  Var_clin = 
  pow(b1CEonSur1Sur2,2) *( pow(bSur1onSur2,2)*sigSqSur2 + SigSqSur1onSur2)+ 
  pow(b2CEonSur1Sur2,2)* sigSqSur2 + SigSqCEonSur1Sur2 + 
  2*b1CEonSur1Sur2*b2CEonSur1Sur2*(bSur1onSur2*sigSqSur2);
  
  R2 = 1 - SigSqCEonSur1Sur2 / Var_clin;
  OptTotal = ((b1CEonSur1Sur2/b2CEonSur1Sur2) * 3 + 3)/12;
}