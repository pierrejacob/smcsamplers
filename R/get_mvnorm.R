#'@rdname get_mvnormal
#'@title Get multivariate normal target
#'@description Define a multivariate normal target
#' with prescribed dimension, mean, covariance.
#'@param dimension dimension of the target
#'@param mean_target mean
#'@param Sigma_target covariance
#'@return list with keys \code{logdensity} and \code{gradlogdensity},
#'@export
get_mvnormal <- function(dimension, mean_target, Sigma_target){
  ### compute precision matrix
  precision <- solve(Sigma_target)
  ### compute Cholesky factor of precision matrix
  precision_chol <- t(chol(precision))
  ### log density of multivariate normal
  logdensity <- function(x){
    return(dmvnorm_cholesky_inverse(x, mean_target, precision_chol))
  }
  ### gradient of log density of multivariate normal
  gradlogdensity <- function(x){
    return(grad_dmvnorm(x, mean_target, precision))
  }

  ### sample from multivariate normal
  rinit <- function(nsamples){
    return(return(rmvn(n = nsamples, mu = mean_target, sigma = Sigma_target)))
  }

  return(list(logdensity = logdensity, gradlogdensity = gradlogdensity, rinit = rinit))
}
