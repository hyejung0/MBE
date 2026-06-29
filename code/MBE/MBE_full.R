#We will produce true treatment effect on the CE in a new RCT 
#with full carry over from previous studies
#Based on MCMC with
#inverse gamma prior on the conditional variances




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


source("../post_theta0_samples_v4.R") #Function 






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







# Full posterior ----------------------------------------------------------



#Do not make diffuse prior on surrogates and keep the intercept
MBE_full_post<-post_theta0_samples_v4(
  mcmc_dat = this.MCMC_dat, 
  sample_dat = sample_dat,
  diffuse=FALSE,
  # diffuse_se = 100,
  # ncores = ncores,
  intercept0=FALSE)


#To save, replace / with __
save_name<-stringr::str_replace_all(study_name,"/","__")
saveRDS(MBE_full_post, file = paste0(save_name, ".rds"))
