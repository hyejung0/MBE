# Appendix: Historical RCTs Model Fit Assessment 
# We show code for how to calculate LOO-CV RMSE , LOO-PD, WAIC, 
#and LOO-CV tolerance interval coverage for the historical model fit. 


rm(list=ls())
library(data.table)



# Read in Data ------------------------------------------------------------

#In ../historical/historical_model_fit.R script,
#we fit the historical model in four different ways.
#We are choosing one arbitrary
#Supposing that its posterior samples were saved as ../historical/historical_fits.rds

myData<-readRDS("../historical/historical_fits.rds")

## LOO-CV RMSE ------------------------------------------------------------

#Leave one out cross validated RMSE
#The MCMC produced draws from the posterior distribution of the held out true clinical effect. 
#To compute the LOO-CV RMSE, these posterior draws were averaged to 
#approximate the posterior mean, then their squared difference from 
#the held out observed clinical endpoint treatment effect was computed. 
#The LOO-CV RMSE then equaled the square root of the average of these squared 
#deviations across studies’ held out clinical effects. 

#Read in the RMSE results
RMSE<-
  mclapply(model_folders,function(model){
    readRDS(paste0("../",model,"/Analysis/LOO-CV_RMSE/LOO_CV_RMSE.rds")) |> list()
  },mc.cores = ncores)
RMSE<-rbindlist(RMSE, idcol = "prior")
setnames(RMSE, "V1","RMSE")
RMSE



## LOO-PD -----------------------------------------------------------------


library(loo)
library(parallel)
#Read in all rds names
LOO_PD<-
  mclapply(model_folders,function(x){
    
    filenames<-c(list.files(paste0("../",x,"/MCMC_output_80000_repar"), pattern = "ALL_STUDIES_loo.rds", full.names = T),
                 list.files(paste0("../",x,"/MCMC_output_80000_c5_repar"), pattern = "ALL_STUDIES_loo.rds", full.names = T)
    )
    #Leave-one-out pointwise log predictive density
    ALL_STUDIES<-readRDS(filenames)
    ALL_STUDIES$estimates["elpd_loo",, drop=FALSE] |> data.table()
    
  }, mc.cores=ncores)
LOO_PD<-rbindlist(LOO_PD, idcol = "prior")
LOO_PD[,LOO_PD:=sprintf("%.1f (SE %.1f)",Estimate, SE)]
LOO_PD


## WAIC -------------------------------------------------------------------

WAIC<-
  mclapply(model_folders,function(x){
    
    filenames<-c(list.files(paste0("../",x,"/MCMC_output_80000_repar"), pattern = "ALL_STUDIES_waic.rds", full.names = T),
                 list.files(paste0("../",x,"/MCMC_output_80000_c5_repar"), pattern = "ALL_STUDIES_waic.rds", full.names = T)
    )
    #Leave-one-oqut pointwise log predictive density
    ALL_STUDIES<-readRDS(filenames)
    ALL_STUDIES$estimates["waic",, drop=FALSE] |> data.table()
    
  }, mc.cores=ncores)
WAIC<-rbindlist(WAIC, idcol = "prior")
WAIC[,WAIC:=sprintf("%.2f (SE %.2f)",Estimate, SE)]
WAIC



## LOO-CV Prl coverage ----------------------------------------------------

#Proportion of cross validation that captured the observed values
capture_observed_values<-
  lapply(model_folders,function(model){
    
    #Read in the folder names that exist in the folder "model"
    #this is a data table showing if each study was captured in the prediction interval or not
    tmp<-readRDS(paste0("../",model,"/Analysis/LOO-PD/capture_observed_value.rds"))
    
    #Calculate the mean captured
    tmp[, .(`95%Prl`=mean(in_95_PI), `90%Prl`=mean(in_90_PI))]
    
  })

coverage<-lapply(capture_observed_values,function(x)apply(x[,.(`95%Prl`,`90%Prl`)],2,mean))
sapply(coverage,"[[", "95%Prl")

PPD_table<-
  data.table(
    prior=names(coverage),
    `95%Prl`=sapply(coverage,"[[", "95%Prl"),
    `90%Prl`=sapply(coverage,"[[", "90%Prl")
  )
PPD_table #Put this in appendix
#             prior    95%Prl    90%Prl
#            <char>     <num>     <num>
# 1:       invgamma 0.9393939 0.8939394
# 2: invgamma_noInt 0.9393939 0.9090909
# 3:       halfnorm 0.9393939 0.8939394
# 4: halfnorm_noInt 0.9393939 0.9090909
(59:65)/66
# 0.8939394 0.9090909 0.9242424 0.9393939 0.9545455 0.9696970 0.9848485
# we see 59/66, 60/66, and 62/66. 


#Check if they are significantly different from nominal 95% and 90% interval 
#using chi-squared test


#For 95\% interval, we should expect 95\% of 66, which is 66*0.95 =62.7 
# we only saw 0.9393, which is 62/66
# So let's test that

alpha<-0.05
N<-66
O<-62
E<-66*(1-alpha)
Obs <- c(O, N-O)
Exp <- c(E, N-E)

# Chi-square without continuity correction
chisq.test(Obs,p=Exp / sum(Exp), rescale.p = TRUE, correct = FALSE)

# Chi-square with continuity correction
chisq.test(Obs,p=Exp / sum(Exp), rescale.p = TRUE, correct = TRUE)

#for both, we get p-value = 0.69 so we can't see a substantial deviation from 0.95.



#For 90\% interval, we should expect 90\% of 66, which is 66*0.90 =59.4 .
#And we got 59 and 60 out of 66
alpha<-0.10
N<-66
E<-66*(1-alpha)
Exp <- c(E, N-E)

# 59
O<-59
Obs <- c(O, N-O)
chisq.test(Obs,p=Exp / sum(Exp), rescale.p = TRUE, correct = FALSE)
chisq.test(Obs,p=Exp / sum(Exp), rescale.p = TRUE, correct = TRUE)
#X-squared = 0.026936, df = 1, p-value = 0.8696
#No problem

# 60
O<-60
Obs <- c(O, N-O)
chisq.test(Obs,p=Exp / sum(Exp), rescale.p = TRUE, correct = FALSE)
chisq.test(Obs,p=Exp / sum(Exp), rescale.p = TRUE, correct = TRUE)
#X-squared = 0.060606, df = 1, p-value = 0.8055

#No problem also!
#No problem!!!!

