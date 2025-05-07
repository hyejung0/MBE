#handcode gibbs sampler that works!


library(data.table)
library(invgamma)



# data --------------------------------------------------------------------


#read in my data
myData<-read.csv("~/Library/CloudStorage/Box-Box/9 Data/0 Latest t99 dataset to use/OriginalModels/nt99a_ce_up_corrs_allacute_20241112.csv")
myData<-data.table(myData)
myData<-myData[ema==1,]
myData<-myData[eventnumber==3,]
myData<-data.table(myData)


anal0 <- read.csv("~/Library/CloudStorage/Box-Box/9 Data/0 Latest t99 dataset to use/OriginalModels/t99_13outcomes_20231201.csv",header=TRUE)       #Full fup/All pts. #there's also data dictionary in the same folder as this. "Data dictionary for txt effects and corrs files.xlsx" file. Actually, "Data dictionary for txt effects and corrs files HPC 20240229.xlsx" file is more comprehensive and has more definition.
anal0 <- anal0[!(is.na(anal0$ev3_loghr1_est)),] #take out rows with missing ev3_loghr1_est values
anal0 <- anal0[!(is.na(anal0$ev3_loghr1_se)),]  #take out rows with missing ev3_loghr1_se
anal0 <- anal0[!(is.na(anal0$gfr_acute_est)),]  #take out rows with missing gfr_acute_est values
anal0 <- anal0[!(anal0$pooled == 2),] #to understand what pooled is, do table(myDataA$pooled_fmt,myDataA$pooled).
anal0<-data.table(anal0)
myData<-
  merge(myData,
        anal0[,.(allrxname,dgfr_chr_corr)],
        all.x=T,
        by="allrxname")
rm(anal0)

myData[,sapply(.SD,mean,na.rm=T),.SDcols = c("beta1_ce_corr","beta3_ce_corr","dgfr_chr_corr")]

#replace the missing values with the mean for beta1_ce_corr,beta3_ce_corr, dgfr_chr_corr
mean_beta1_ce_corr<-myData[,mean(beta1_ce_corr, na.rm=T)] #We impute with the mean
mean_beta3_ce_corr<-myData[,mean(beta3_ce_corr, na.rm=T)]
mean_dgfr_chr_corr<-myData[,mean(dgfr_chr_corr, na.rm=T)]

#Keet the originals first 
myData[,Orig_beta1_ce_corr:=beta1_ce_corr]
myData[,Orig_beta3_ce_corr:=beta3_ce_corr]
myData[,Orig_dgfr_chr_corr:=dgfr_chr_corr]

#Impute
myData[is.na(beta1_ce_corr),beta1_ce_corr:=mean_beta1_ce_corr]
myData[is.na(beta3_ce_corr),beta3_ce_corr:=mean_beta3_ce_corr]
myData[is.na(dgfr_chr_corr),dgfr_chr_corr:=mean_dgfr_chr_corr]



# "Sur2Est", "Sur1Est", "ClnEst","ClnSE","Sur1SE","Sur2SE","R1Clin","R2Clin","R12"
dat_list<-data.table(
  Sur1Est = myData$beta31_beta32_est,
  Sur1SE  = myData$beta31_beta32_se,

  Sur2Est = myData$zbeta11_beta12_est,
  Sur2SE  = myData$zbeta11_beta12_se,

  ClnEst  = myData$loghr1_est,
  ClnSE   = myData$loghr1_se,

  R1Clin = myData$beta3_ce_corr,
  R2Clin = myData$beta1_ce_corr,
  R1R2   = myData$dgfr_chr_corr
)








# initial staring values --------------------------------------------------


#Define initial values

init_list<-list(
  #this is same as Tom's SAS code
  list(
    "b2CEonSur1Sur2"=0,
    "SigSqCEonSur1Sur2"=0.000010,
    "alphaCEonSur1Sur2"=0,
    "b1CEonSur1Sur2"=0,
    
    "SigSqSur1onSur2"=1,
    "alphaSur1onSur2"=0,
    "bSur1onSur2"=0,
    
    "muSur2"=-.4645,
    "sigSqSur2"=2.1993
  ),
  
  
  #this is at the lower value
  list(
    "alphaCEonSur1Sur2"=-100, #One standard deviation below
    "b2CEonSur1Sur2"=-100, #One standard deviation below
    "b1CEonSur1Sur2"=-100, #One standard deviation below
    "SigSqCEonSur1Sur2"=qinvgamma(0.25, shape=0.261,rate=0.000408), #this is the cut off value for low heterogeneity for SD
    
    "alphaSur1onSur2"=-100, #One standard deviation below
    "bSur1onSur2"=-100, #One standard deviation below
    "SigSqSur1onSur2"=qinvgamma(0.25, shape=0.261,rate=0.005),#this is the cut off value for low heterogeneity for SD
    
    "muSur2"=-100, #One standard deviation below
    "sigSqSur2"=qinvgamma(0.25, shape=0.261,rate=0.000408)#this is the cut off value for low heterogeneity
  ),
  
  
  #this is at the middle value
  list(
    "alphaCEonSur1Sur2"=0, #mean
    "b2CEonSur1Sur2"=0, #mean
    "b1CEonSur1Sur2"=0, #mean
    "SigSqCEonSur1Sur2"=qinvgamma(0.5, shape=0.261,rate=0.000408), #median
    
    "alphaSur1onSur2"=0,
    "bSur1onSur2"=0,
    "SigSqSur1onSur2"=qinvgamma(0.5, shape=0.261,rate=0.005),#this is the cut off value for low heterogeneity for SD
    
    "muSur2"=0, #One standard deviation below
    "sigSqSur2"=qinvgamma(0.5, shape=0.261,rate=0.000408)#this is the cut off value for low heterogeneity
  ),
  
  
  #this is at the higher value
  list(
    "alphaCEonSur1Sur2"= 100, #One standard deviation above
    "b2CEonSur1Sur2"= 100, #One standard deviation above
    "b1CEonSur1Sur2"= 100, #One standard deviation above
    "SigSqCEonSur1Sur2"=qinvgamma(0.75, shape=0.261,rate=0.000408), #this is the cut off value for low heterogeneity for SD
    
    "alphaSur1onSur2"= 100, #One standard deviation above
    "bSur1onSur2"= 100, #One standard deviation above
    "SigSqSur1onSur2"=qinvgamma(0.75, shape=0.261,rate=0.005),#this is the cut off value for medium heterogeneity for SD
    
    "muSur2"= 100, #One standard deviation above
    "sigSqSur2"=qinvgamma(0.75, shape=0.261,rate=0.000408) 
  )
  
)



# Gibbs sampling ----------------------------------------------------------


source("gibbs_sampler.R") #this one uses conditional postierior distribution individually. 
debugonce(gibbs_sampler)
set.seed(30)
gibbs_out<-gibbs_sampler(dat_list,inits = init_list,burn.out = 0.5,n_iter = 30000,parallel = F,ncores = 5)
gibbs_out$quantile.par
gibbs_out$Rhat
gibbs_out$effective.size

source("check.mixing.R")
check.mixing(as.matrix(gibbs_out$draws), nChains = 4)
library(bayesplot)

#For each parameters, make them into an array of matrix where column is chain.

dt<-copy(gibbs_out$draws)
dt[, iteration := seq_len(.N), by = chains]

# dims:
n_iter  <- max(dt$iteration)
chains  <- sort(unique(dt$chains))
n_chain <- length(chains)
pars    <- param_cols
n_par   <- length(pars)

# build an empty array
arr <- array(
  NA_real_,
  dim = c(n_iter, n_chain, n_par),
  dimnames = list(
    iteration = 1:n_iter,
    chain     = as.character(chains),
    parameter = pars
  )
)

# fill it
for (j in seq_along(chains)) {
  ch <- chains[j]
  for (k in seq_along(pars)) {
    arr[, j, k] <- dt[chains == ch, get(pars[k])]
  }
}

# trace plot
pdf("gibbs_individiual_traceplot.pdf",width=15,height=10)
mcmc_trace(arr)
dev.off()
#I think pdf is inch


# Rstan -------------------------------------------------------------------

library(data.table)
library(rstan)
library(tictoc)
library(rtf)
library(bayesplot)
library(ggplot2)


## Data -------------------------------------------------------------------



# With multi-arm corrs:


data_sub = list(
  
  nStudies=dim(myData)[1],
  
  chrn_eff=myData$beta31_beta32_est,
  chrn_se=myData$beta31_beta32_se,
  
  actslp_eff=myData$zbeta11_beta12_est,
  actslp_se=myData$zbeta11_beta12_se,
  
  # acr_eff=acr_eff,
  # acr_se=acr_se,
  
  ClnEst=myData$loghr1_est,
  ClnSE=myData$loghr1_se,
  
  Rclnchron=myData$beta3_ce_corr,
  Rclnacute=myData$beta1_ce_corr,
  # Rclnacr=Rclnacr,
  Rchronacute=myData$dgfr_chr_corr
  # Rchronacr=Rchronacr,
  # Racuteacr=Racuteacr
  
)  

data_sub_dt<-as.data.table(data_sub)

#Set up initial values same as ones used in Gibbs sampling
init_list

## Setup RStan parameters -------------------------------------------------



options(mc.cores = min(length(init_list),5))
num.iter=30000
num.warmup=num.iter/2
num.chains<-length(init_list)
num.thin=1









## Define variables of interests ------------------------------------------


all_pars<-
  c(
    "b2CEonSur1Sur2",
    "SigSqCEonSur1Sur2",
    "alphaCEonSur1Sur2",
    "b1CEonSur1Sur2",
    
    "SigSqSur1onSur2",
    "alphaSur1onSur2",
    "bSur1onSur2",
    
    "muSur2",
    "sigSqSur2",
    "R2",
    "OptTotal"
  )


set.seed(1)
tic()
myStanmod = stan("chr_acute_hyejung.stan",
                 data=data_sub,
                 iter=num.iter,
                 warmup=num.warmup,
                 chains=num.chains,
                 thin=num.thin,
                 init = init_list,
                 control=list(adapt_delta=0.99,max_treedepth=25))

toc()

#save posterior samples
posterior_samples<-extract(myStanmod,all_pars,permuted=TRUE)
posterior_samples<-lapply(posterior_samples,unlist)
posterior_samples<-do.call(cbind,posterior_samples)
posterior_samples<-data.table(posterior_samples)

quantile_par<-apply(posterior_samples,2,function(x)quantile(x,c(.025,.25,.5,.75,.975)))
quantile_par<-t(quantile_par)
quantile_par




# Compare Gibbs vs RStan --------------------------------------------------

#Rstan
quantile_par[all_pars,]

#Gibbs
gibbs_out$quantile.par[all_pars,]

compare_two<-
lapply(all_pars,function(par){
  
  tmp<-rbind(
    "rstan"=quantile_par[par,, drop=F],
    "Gibbs"=gibbs_out$quantile.par[par,,drop=F]
  )
  rownames(tmp)<-c("rstan","Gibbs")
  
  tmp
})
names(compare_two)<-all_pars
compare_two


#very different for 
#b2CEonSur1Sur2
#alphaCEonSur1Sur2
#alphaSur1onSur2
#bSur1onSur2
#sigSqSur2



#Slightly different for:
#b1CEonSur1Sur2
#SigSqSur1onSur2
#OptTotal



# Gibbs multisampling -----------------------------------------------------

source("gibbs_sampler_v2.R")
# debugonce(gibbs_sampler_v2)
set.seed(30)
tic()
gibbs_out2<-gibbs_sampler_v2(dat_list,inits = init_list,burn.out = 0.5,n_iter = 70000,parallel = T,ncores = 5)
toc()
#104.976 sec elapsed

dt<-copy(gibbs_out2$draws)
dt[, iteration := seq_len(.N), by = chains]

# dims:
n_iter  <- max(dt$iteration)
chains  <- sort(unique(dt$chains))
n_chain <- length(chains)
pars    <- all_pars
n_par   <- length(pars)

# build an empty array
arr <- array(
  NA_real_,
  dim = c(n_iter, n_chain, n_par),
  dimnames = list(
    iteration = 1:n_iter,
    chain     = as.character(chains),
    parameter = pars
  )
)

# fill it
for (j in seq_along(chains)) {
  ch <- chains[j]
  for (k in seq_along(pars)) {
    arr[, j, k] <- dt[chains == ch, get(pars[k])]
  }
}

# trace plot
pdf("gibbs_MH_traceplot.pdf",width=15,height=10)
mcmc_trace(arr)
dev.off()


gibbs_out2$quantile.par
gibbs_out2$Rhat
gibbs_out2$effective.size

compare_two2<-
  lapply(all_pars,function(par){
    
    tmp<-rbind(
      "rstan"=quantile_par[par,, drop=F],
      "Gibbs"=gibbs_out2$quantile.par[par,,drop=F]
    )
    rownames(tmp)<-c("rstan","Gibbs")
    
    tmp
  })
names(compare_two2)<-all_pars

compare_two2


# $b2CEonSur1Sur2
#              2.5%         25%         50%         75%       97.5%
# rstan -0.03493016 -0.02961140 -0.02693208 -0.02424201 -0.01905605
# Gibbs -0.03498580 -0.02962127 -0.02690092 -0.02425828 -0.01907274
# 
# $SigSqCEonSur1Sur2
#               2.5%         25%         50%         75%      97.5%
# rstan 0.0002390212 0.001360081 0.003439827 0.007051433 0.01791135
# Gibbs 0.0002355252 0.001326342 0.003448320 0.006999421 0.01786817
# 
# $alphaCEonSur1Sur2
#              2.5%         25%         50%          75%      97.5%
# rstan -0.09679219 -0.05308773 -0.03142482 -0.010637486 0.03095102
# Gibbs -0.09983625 -0.05365410 -0.03144087 -0.009979343 0.03334812
# 
# $b1CEonSur1Sur2
#             2.5%        25%        50%        75%      97.5%
# rstan -0.4299655 -0.3719212 -0.3442317 -0.3167331 -0.2639762
# Gibbs -0.4304283 -0.3724593 -0.3441372 -0.3168041 -0.2631232
# 
# $SigSqSur1onSur2
#            2.5%       25%       50%       75%     97.5%
# rstan 0.2404655 0.3276211 0.3864989 0.4584370 0.6407811
# Gibbs 0.2398827 0.3274984 0.3867748 0.4579864 0.6453269
# 
# $alphaSur1onSur2
#            2.5%       25%       50%       75%     97.5%
# rstan 0.3313320 0.4588773 0.5244101 0.5903673 0.7206615
# Gibbs 0.3276721 0.4571600 0.5237828 0.5903031 0.7216720
# 
# $bSur1onSur2
#              2.5%         25%          50%         75%      97.5%
# rstan -0.03982998 -0.01852746 -0.007643443 0.003093673 0.02472885
# Gibbs -0.03994162 -0.01864927 -0.007681486 0.003261213 0.02488523
# 
# $muSur2
#            2.5%       25%       50%        75%     97.5%
# rstan -3.041431 -1.930883 -1.356439 -0.7829743 0.3726744
# Gibbs -3.048646 -1.931939 -1.358147 -0.7791727 0.3555124
# 
# $sigSqSur2
#           2.5%      25%      50%      75%    97.5%
# rstan 25.39144 33.30261 38.63897 45.03615 60.86378
# Gibbs 25.36164 33.29268 38.61034 44.94033 60.87276
# 
# $R2
#            2.5%      25%       50%       75%     97.5%
# rstan 0.7763689 0.908017 0.9537871 0.9811801 0.9967556
# Gibbs 0.7765683 0.908601 0.9535645 0.9816968 0.9969086
# 
# $OptTotal
#           2.5%      25%      50%      75%    97.5%
# rstan 2.605155 3.133002 3.443984 3.813924 4.779707
# Gibbs 2.594431 3.130913 3.446326 3.812032 4.774016




# Update ------------------------------------------------------------------

# now we want to do some update 
# Suppose we observe a new data.

