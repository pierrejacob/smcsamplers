#'@rdname get_student_t
#'@title Get student t distribution
#'@description Define a student t distribution
#' with prescribed dimension, mean, variance.
#'@param dimension dimension of the target
#'@param degree degrees of freedom
#'@param mean_target mean
#'@param Sigma_target covariance
#'@return list with keys \code{logdensity} and \code{gradlogdensity},
#'@export
get_student_t <- function(dimension, degree, mean_proposal, Sigma_proposal){
  ### compute precision matrix
  precision <- solve(Sigma_proposal)

  ### compute Cholesky factor of precision matrix
  precision_chol <- t(chol(precision))

  ### log density of multivariate student t distribution
  factor1 <- 0.5 * (dimension + degree)
  constant <- lgamma(factor1) - lgamma(0.5 * degree) -
              0.5 * dimension * ( log(degree) + log(pi) ) + sum(log(diag(precision_chol)))
  logdensity <- function(x){
    return(dmvstudent_t_cholesky_inverse(x, degree, mean_proposal, precision_chol, constant, factor1))
  }

  ### gradient of log density of multivariate student t
  factor2 <- - (degree + dimension) / degree
  gradlogdensity <- function(x){
    return(grad_dmvstudent_t(x, degree, mean_proposal, precision, factor2))
  }

  ### sample from multivariate student t
  rinit <- function(nsamples){
    return(rmt(nsamples, mean = mean_proposal, S = Sigma_proposal, df = degree))
  }

  return(list(logdensity = logdensity, gradlogdensity = gradlogdensity, rinit = rinit))
}
