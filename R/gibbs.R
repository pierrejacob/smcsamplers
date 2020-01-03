#' @rdname gibbs
#' @title Run the Metropolis-within-Gibbs algorithm
#' @param target list with keys: \code{logdensity} evaluates log target density, \code{gradlogdensity} returns its gradient
#' @param proposal list with keys: \code{rinit} samples from proposal, \code{logdensity} evaluates log proposal density, \code{gradlogdensity} returns its gradient
#' @param gibbs_parameters a vector of step sizes used in the Metroplis algorithm
#' @param gamma_current the current inverse temperature
#' @export

gibbs <- function(target, proposal, gibbs_parameters, gamma_current){
  ### a vector of step sizes used in the Metropolis algorithm
  stepsize <- gibbs_parameters$stepsize

  ### construct the Metropolis-within-Gibbs kernel
  kernel <- function(chain_state, logtarget_current, logproposal_current){
    dimension <- ncol(chain_state)
    nparticles <- nrow(chain_state)
    accept_ratio <- matrix(NA, nrow = nparticles, ncol = dimension)
    ### Gibbs sampling with systematic scan
    for(j in 1:dimension){
      forward <- chain_state
      ### propose a new jth coordinate
      forward[, j] <- forward[, j] + rnorm(nparticles, mean = 0, sd = stepsize[j])
      logtarget_forward <- target$logdensity(forward)
      logproposal_forward <- proposal$logdensity(forward)
      loglik_forward <- gamma_current * logtarget_forward + (1 - gamma_current) * logproposal_forward
      loglik_current <- gamma_current * logtarget_current + (1 - gamma_current) * logproposal_current
      ### calculate the log acceptance ratio
      accept_ratio[, j] <- loglik_forward - loglik_current
      ### the Metropolis acceptance or rejection step
      accept_logical <- log(runif(nparticles)) < accept_ratio[, j]
      chain_state[accept_logical, ] <- forward[accept_logical, ]
      logtarget_current[accept_logical] <- logtarget_forward[accept_logical]
      logproposal_current[accept_logical] <- logproposal_forward[accept_logical]
    }

    return(list(chain_state = chain_state, logtarget = logtarget_current, logproposal = logproposal_current,
                accept_ratio = mean(pmin(1, exp(accept_ratio)))))
  }
  return(list(kernel = kernel))
}
