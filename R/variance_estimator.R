## From the output of SMCS, compute asymptotic variance estimator for the normalizing constant,
## as in Lee & Whiteley 2018 and as described in Section 4.1 of the article.
#'@export
variance_estimator <- function(results){
  nparticles <- results$particles$n
  roots <- results$roots_history[[length(results$roots_history)]]
  np1 <- length(results$lambdas)
  prod_ <- (nparticles/(nparticles-1))^(np1)
  count_sameancestor <- rep(0, nparticles)
  for (i in 1:nparticles){
    count_sameancestor[roots[i]] <- count_sameancestor[roots[i]] + 1
  }
  var_estim <- (1 - prod_) + prod_ * sum((count_sameancestor/nparticles)^2)
  return(var_estim)
}
