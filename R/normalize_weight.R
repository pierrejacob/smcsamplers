#'@rdname normalize_weight
#'@title Normalize vector of log-weights
#'@description Takes a vector of real values, and return the vector of normalized exponentiated values,
#' as well as average of the exponentiated unnormalized values.
#'@param logweights is a vector of n real values
#'@return a vector of n non-negative values summing to one
#'@examples
#' N <- 1000
#' logweights <- rnorm(N)
#' normalize_weight_results <- normalize_weight(logweights)
#' normalized_weights <- normalize_weight_results$nw
#'@export
normalize_weight <- function(logweights){
  mlw <- max(logweights)
  avew <- mlw + log(mean(exp(logweights - mlw))) # weight average
  w <- exp(logweights - mlw)
  nw <- w / sum(w) # normalized weights
  return(list(nw = nw, avew = avew))
}
