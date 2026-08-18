#' Matrix Mean and Variance
#'
#' @description Generate mean and covariance using prior data. We will use this to construct weights to update prior \eqn{\xi}.
#'
#' @param MCMC_dat A data.table object containing MCMC samples.
#' @param diffuse_se A numeric value specifying the standard error for the diffuse prior.
#' @param diffuse A logical value indicating whether to use a diffuse prior. Defaults to `TRUE`.
#' @param intercept0 A logical value indicating whether to center the intercept term of the meta-regression to zero. Defaults to `TRUE`.
#'
#' @return A list containing:
#' \describe{
#'    \item{means}{A list of mean vectors for each MCMC sample.}
#'    \item{covar}{A list of covariance matrices for each MCMC sample.}
#'    \item{dt}{A data.table containing the mean and covariance values for each MCMC sample.}
#' }
#'
#' @examples
#' # Use the example MCMC data provided in the package `historical_posterior`
#' result_matrix<-vector_to_matrix(
#' MCMC_dat = historical_posterior,
#' diffuse_se = 100,
#' diffuse = TRUE,
#' intercept0 = TRUE
#' )
#'
#' head(result_matrix)
#' @noRd
vector_to_matrix<-function(MCMC_dat,diffuse_se=100, diffuse=TRUE, intercept0=TRUE){


  #Make sure we are using data.table format
  this.MCMC_dat<-data.table::data.table(MCMC_dat)

  #Save number of MCMC samples
  B<-nrow(this.MCMC_dat)


  #If we are specifying diffuse prior, set the corresponding MCMC parameters to what we want to set.
  if(diffuse){
    this.MCMC_dat[,bSur1onSur2:=0]
    this.MCMC_dat[,SigSqSur1onSur2:=diffuse_se^2]
    this.MCMC_dat[,sigSqSur2:=diffuse_se^2]
    this.MCMC_dat[,muSur2:=0]
    this.MCMC_dat[,alphaSur1onSur2:=0]
  }

  #If we want to shift the mean of the intercept to 0
  if(intercept0){
    this.MCMC_dat[,alphaCEonSur1Sur2:=scale(alphaCEonSur1Sur2, center=TRUE, scale=FALSE)]
  }

  #Variance
  #define variance structure
  this.MCMC_dat[,Sigma33:=sigSqSur2]
  this.MCMC_dat[,Sigma22:=bSur1onSur2^2 * sigSqSur2 + SigSqSur1onSur2]
  this.MCMC_dat[,Sigma23:=bSur1onSur2*sigSqSur2]
  this.MCMC_dat[,Sigma13:=b1CEonSur1Sur2*(bSur1onSur2*sigSqSur2) + b2CEonSur1Sur2*sigSqSur2]
  this.MCMC_dat[,Sigma12:=b1CEonSur1Sur2*(bSur1onSur2^2*sigSqSur2 + SigSqSur1onSur2)+ b2CEonSur1Sur2*bSur1onSur2*sigSqSur2]
  this.MCMC_dat[,Sigma11:=b1CEonSur1Sur2^2 *(bSur1onSur2^2*sigSqSur2 + SigSqSur1onSur2)+b2CEonSur1Sur2^2 * sigSqSur2 + SigSqCEonSur1Sur2 + 2*b1CEonSur1Sur2*b2CEonSur1Sur2*(bSur1onSur2*sigSqSur2)]



  #Array => list is fastest without using parallel



  #An array of covariances
  these_cols <- c("Sigma11","Sigma12","Sigma22","Sigma13","Sigma23","Sigma33")
  tmp <- as.matrix(this.MCMC_dat[, ..these_cols])   # B x 6
  covars_array <- array(0, dim = c(3, 3, B))
  # Fill upper triangle incl diag
  ut <- upper.tri(matrix(0, 3, 3), diag = TRUE)
  for (i in seq_len(B)) {
    covars_array[,,i][ut] <- tmp[i, ]
    covars_array[,,i][!ut] <- t(covars_array[,,i])[!ut]
  }
  #Convert to a list
  covars <- lapply(seq_len(dim(covars_array)[3]), function(i) covars_array[,,i])

  #Do likewise for means
  this.MCMC_dat[,mu1:=alphaCEonSur1Sur2 + b1CEonSur1Sur2*(alphaSur1onSur2 + bSur1onSur2*muSur2) + b2CEonSur1Sur2*muSur2]
  this.MCMC_dat[,mu2:=alphaSur1onSur2 + bSur1onSur2*muSur2]
  this.MCMC_dat[,mu3:=muSur2]
  tmp<-as.matrix(this.MCMC_dat[,.(mu1, mu2, mu3)])
  means_array <- array(0, dim = c(3, 1, B))
  for (i in seq_len(B)) {
    means_array[,,i] <- tmp[i, ]
  }
  means <- lapply(seq_len(dim(means_array)[3]), function(i) means_array[,,i,drop=FALSE])


  return(list(mean=means, covar=covars, dt=this.MCMC_dat[,.SD,.SDcols=c("mu1","mu2","mu3",these_cols)]))


}


