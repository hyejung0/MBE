

# MCMC data set -----------------------------------------------------------

# 1. Read the .rds file into a named R object
historical_posterior <- readRDS("data-raw/ALL_STUDIES.rds")

#Too many datasets. I will only save the last 1,000 samples from each chain.
historical_posterior<-data.table(historical_posterior)
max_iter<-historical_posterior[,max(.iteration)]
all_pars<-c(
  "alphaCEonSur1Sur2" ,
  "b1CEonSur1Sur2",
  "b2CEonSur1Sur2",
  "SigSqCEonSur1Sur2",

  "alphaSur1onSur2",
  "bSur1onSur2",
  "SigSqSur1onSur2",

  "muSur2",
  "sigSqSur2"
)
historical_posterior<-historical_posterior[.iteration>(max_iter - 1000), ..all_pars]

# 2. Save it as an .rda file into the official data/ folder
usethis::use_data(historical_posterior, overwrite = TRUE)

# # #Generates bzip2 file because the raw .rda file is too big
# usethis::use_data(historical_posterior,compress ="bzip2", overwrite = TRUE)

# simulated dataset -------------------------------------------------------

#Generated 66 example trials from Yizhen's code.
trial_sim_dat<-read.csv("data-raw/MBE_simulated_66_trials.csv")
trial_sim_dat<-data.table::data.table(trial_sim_dat)
# Save it as an .rda file into the official data/ folder
usethis::use_data(trial_sim_dat, overwrite = TRUE)

# simulated dataset for interim analysis ----------------------------------

# 1. Read the .rds file into a named R object
interim_sim_dat <- readRDS("data-raw/dat.rds")

# Keep only the single clinical endpoint and the two retained surrogate
# endpoints (chronic and acute eGFR slopes). Third-surrogate fields are
# intentionally excluded from both the raw and packaged example data.
interim_sim_dat <- interim_sim_dat[, c(
  "case", "analysis.month",
  "R1Clin", "R2Clin", "R12",
  "ClnEst", "ClnSE",
  "Sur2Est", "Sur2SE",
  "Sur1Est", "Sur1SE"
), with = FALSE]

# 2. Save it as an .rda file into the official data/ folder
usethis::use_data(interim_sim_dat, overwrite = TRUE)

