### script that defines some target distributions: 'banana' and multivariate Normals

## banana distribution
#'@export
get_banana <- function(){
  # target log-density evaluation
  bananatarget <- function(x) -(1-x[1])^2 - 10*((x[2]-x[1]^2)^2)
  # gradient of target log-density
  bananagradtarget <- function(x) c(-2*(x[1]-1) + 40*x[1]*(x[2]-x[1]^2),
                                    -20 * (x[2]-x[1]^2))
  # function to compute gradient and log pdf for each column of 'xs'
  bananaeval <- function(xs, ...){
    return(list(gradlogpdf = apply(xs, 2, bananagradtarget), logpdf = apply(xs, 2, bananatarget)))
  }
  return(list(eval = bananaeval))
}

## multivariate Normal distribution
#'@export
get_mvnormal <- function(mean, variance){
  dimension <- length(mean)
  precision <- solve(variance)
  cholvariance <- chol(variance)
  cholinvvariance <- t(chol(solve(variance)))
  eval_ <- function(xs, ... ){
    return(list(gradlogpdf = t(smcsamplers:::grad_dmvnorm_precision(t(xs), mean = mean, precision)),
                logpdf = smcsamplers:::dmvnorm_cholesky_inverse(t(xs), mean, cholinvvariance)))

  }
  return(list(eval = eval_,
              generate = function(n) t(smcsamplers:::rmvnorm_cholesky_(n, mean, cholvariance))))
}

## multivariate Normal with diagonal covariance matrix (faster than 'get_mvnormal')
#'@export
get_mvnormal_diag <- function(mean, variance){
  sd_ <- sqrt(variance)
  eval_ <- function(xs, ... ){
    return(list(gradlogpdf = -(xs - mean) / variance,
                logpdf = colSums(dnorm((xs - mean), mean = 0, sd = sd_, log = TRUE))))
  }
  return(eval_)
}

