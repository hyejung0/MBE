# This file is a R version of the "Frequentist 2 Surrogates 2020-Hyejung.sas" written by Tom
# Conversion was done by ChatGPT.
# This script demonstrates estimateing MLE for parameters using all 66 

rm(list=ls())
library(data.table)     
library(mvmeta)    # for multivariate meta‐analysis


# Data --------------------------------------------------------------------


myData <- readRDS("evt3data.rds")
myData<-data.table(myData)

#Keep only the relevant variable
myData<-
myData[
  ,
  .(
    allrx,
    beta31_beta32_est,  beta31_beta32_se,   # y2  & SE for chronic‐slope surrogate
    zbeta11_beta12_est, zbeta11_beta12_se,  # y3  & SE for acute‐slope surrogate
    loghr1_est,         loghr1_se,          # y1  & SE for clinical endpoint (log‐HR)
    beta3_ce_corr,      beta1_ce_corr,      # R1Clin (corr: clinical vs chronic), R2Clin (corr: clinical vs acute)
    dgfr_chr_corr       # R1R2 (corr: chronic vs acute)
  )
]

#Rename the variables
setnames(
  myData,
  old=c(
    "beta31_beta32_est",
    "beta31_beta32_se",
    "zbeta11_beta12_est",
    "zbeta11_beta12_se",
    "loghr1_est",
    "loghr1_se",
    "beta3_ce_corr",    
    "beta1_ce_corr",    
    "dgfr_chr_corr"     
  ),
  new=c(
    "y2",   
    "y2_se",
    "y3",   
    "y3_se",
    "y1",   
    "y1_se",
    "r12",  # corr(y1,y2)
    "r13",  # corr(y1,y3)
    "r23"   # corr(y2,y3)
  )
)



# ──────────────────────────────────────────────────────────────────────────────
# 4. for each study (allrx), compute the 3×3 sampling‐covariance matrix C_i:
#      C_i = [ c11  c12  c13 ;
#              c12  c22  c23 ;
#              c13  c23  c33 ]
#    where c11 = SE(y1)^2, c22 = SE(y2)^2, c33 = SE(y3)^2,
#          c12 = r12 * sqrt(c11 * c22), etc.
# ──────────────────────────────────────────────────────────────────────────────

myData[,c11 := (y1_se)^2]
myData[,c22 := (y2_se)^2]
myData[,c33 := (y3_se)^2]
myData[ ,c12 := r12 * sqrt(c11 * c22)]
myData[ ,c13 := r13 * sqrt(c11 * c33)]
myData[ ,c23 := r23 * sqrt(c22 * c33)]




# ──────────────────────────────────────────────────────────────────────────────
# 5. build a list of sampling‐covariance matrices, one 3×3 per study
#    (mvmeta expects “S_list” to be a list-of‐matrices)
# ──────────────────────────────────────────────────────────────────────────────
# We’ll split 'myData' by allrx; each element is a single‐row data.table.


#An array of covariances
these_cols <- c("c11","c12","c22","c13","c23","c33")
tmp <- as.matrix(myData[, ..these_cols])   # nrow(myData) x 6
covars_array <- array(0, dim = c(3, 3, nrow(myData)))

# Fill upper triangle incl diag
ut <- upper.tri(matrix(0, 3, 3), diag = TRUE)
for (i in seq_len(nrow(myData))) {
  covars_array[,,i][ut] <- tmp[i, ]
  covars_array[,,i][!ut] <- t(covars_array[,,i])[!ut]
}
#Convert to a list
S_list <- lapply(seq_len(dim(covars_array)[3]), function(i) covars_array[,,i])



# ──────────────────────────────────────────────────────────────────────────────
# 6. build a matrix Y of observed effect‐estimates (y1,y2,y3) per study
#    mvmeta wants a k×p matrix, where k = #studies, p = #outcomes (here p=3)
# ──────────────────────────────────────────────────────────────────────────────
Y <- myData[,.(y1,y2,y3)] |> as.matrix()


# ──────────────────────────────────────────────────────────────────────────────
# 7. fit the ML‐based random‐effects model exactly like PROC NLMIXED did:
#      u_i = (u1_i, u2_i, u3_i)' ~ MVN( (mu1,mu2,mu3)', PSI )
#      y_i | u_i ~ MVN( u_i, C_i )
#    (i=1,…,N studies).  mvmeta() will estimate:
#      • the vector of means (mu1,mu2,mu3)
#      • the between‐study covariance matrix PSI (3×3)
#
#    “method='ml'” matches SAS’s MLE.  By default mvmeta uses X=Intercept only,
#    so we get exactly the same structure as a “no‐covariate” random‐effects model.
# ──────────────────────────────────────────────────────────────────────────────
fit <- mvmeta(
  Y ~ 1,
  S = S_list,
  method = "ml"
)

summary(fit)


# ──────────────────────────────────────────────────────────────────────────────
# 8. extract the estimated between‐study covariance (PSI) and the means
#    (these correspond to the PROC NLMIXED “parms mu1=… mu2=… mu3=…” output)
# ──────────────────────────────────────────────────────────────────────────────
mu_hat  <- coef(fit)      # a length‐3 vector c(mu1, mu2, mu3)
PSI_hat <- fit$Psi        # 3×3 between‐study covariance matrix

# check:
#   mu_hat[1] ≈ mu1 •
#   mu_hat[2] ≈ mu2 •
#   mu_hat[3] ≈ mu3

# Similarly, PSI_hat[1,1] ≈ psi11, PSI_hat[2,2] ≈ psi22, etc.


# ──────────────────────────────────────────────────────────────────────────────
# 9. compute the regression “true effects” of u1 on (u2,u3):
#    SAS used:
#      DetPSI = psi22*psi33 − psi23^2
#      β2 = (ψ12·ψ33 − ψ13·ψ23)/DetPSI
#      β3 = (ψ13·ψ22 − ψ12·ψ23)/DetPSI
#    which are the coefficients from
#      u1_i | (u2_i, u3_i) ~ Normal( μ1 + β2 (u2_i−μ2) + β3 (u3_i−μ3),  variance=… )
#
#    Here’s the same in R:
# ──────────────────────────────────────────────────────────────────────────────
psi11 <- PSI_hat[1, 1]
psi22 <- PSI_hat[2, 2]
psi33 <- PSI_hat[3, 3]
psi12 <- PSI_hat[1, 2]
psi13 <- PSI_hat[1, 3]
psi23 <- PSI_hat[2, 3]

DetPSI <- psi22 * psi33 - psi23^2

ipsi22 <- psi33/DetPSI
ipsi33 <- psi22/DetPSI
ipsi23 <- -1*psi23/DetPSI

beta2 <- (psi12 * psi33 - psi13 * psi23) / DetPSI
beta3 <- (psi13 * psi22 - psi12 * psi23) / DetPSI

# Proportion of true‐variance explained (R²_true) in SAS:
#   pvarTru = ψ11 − (ψ12·β2 + ψ13·β3)
#   pvarTruAlt = β2·ψ22·β2 + 2·β2·ψ23·β3 + β3·ψ33·β3
pvarTru     <-(psi12*ipsi22*psi12 + psi12*ipsi23*psi13 + psi13*ipsi23*psi12 + psi13*ipsi33*psi13)
pvarTruAlt  <- beta2*psi22*beta2 +
  2*beta2*psi23*beta3 +
  beta3*psi33*beta3

R2_true     <- pvarTru    / psi11
R2_true_alt <- pvarTruAlt / psi11
RMSE_tru    <- sqrt(psi11 - pvarTru)

# ──────────────────────────────────────────────────────────────────────────────
# 10. print out exactly the same “estimates” that SAS PROC NLMIXED did:
#     • intercept for the regression of u1 on (u2,u3): μ1 − μ2·β2 − μ3·β3
#     • “Beta Sur 1” (= β2), “Beta Sur 2” (= β3)
#     • “R2” and “R2 Alt” and “RMSETru”
# ──────────────────────────────────────────────────────────────────────────────
Intercept_true <- mu_hat[1] - mu_hat[2]*beta2 - mu_hat[3]*beta3

# compile into a small table
results <- data.frame(
  Parameter       = c("Intercept (μ1 − μ2·β2 − μ3·β3)",
                      "Beta Surrogate 1 (β2)",
                      "Beta Surrogate 2 (β3)",
                      "R2 (True Proportion)",
                      "R2 Alt",
                      "RMSE_True"),
  Estimate        = c(Intercept_true,
                      beta2,
                      beta3,
                      R2_true,
                      R2_true_alt,
                      RMSE_tru)
)

print(results, row.names = FALSE)
#                      Parameter    Estimate
# Intercept (μ1 − μ2·β2 − μ3·β3) -0.03522787
#          Beta Surrogate 1 (β2) -0.34302303
#          Beta Surrogate 2 (β3) -0.02693480
#           R2 (True Proportion)  0.93682812
#                         R2 Alt  0.93682812
#                      RMSE_True  0.06702862

