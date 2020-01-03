#' @rdname slice_gibbs
#' @title Run slice Gibbs sampler
#' @param target list with keys: \code{logdensity} evaluates log target density, \code{gradlogdensity} returns its gradient
#' @param proposal list with keys: \code{rinit} samples from proposal, \code{logdensity} evaluates log proposal density, \code{gradlogdensity} returns its gradient
#' @param gamma_current the current inverse temperature
#' @param conditional_distribution the conditional sampling distribution
#' @export

slice_gibbs <- function(target, proposal, conditional_distribution, gamma_current){
  cond_x <- conditional_distribution$conditional_x
  cond_u <- conditional_distribution$conditional_u

  kernel <- function(chain_state){
    ### chain_state: N by p matrix
    dimension <- ncol(chain_state)
    nparticles <- nrow(chain_state)

    ### a systematic scan
    for(i in 1:dimension){
      u <- cond_u(chain_state, i, gamma_current)
      chain_state[, i] <- cond_x(chain_state, u, i, gamma_current)
    }

    logtarget <- target$logdensity(chain_state)
    logproposal <- proposal$logdensity(chain_state)

    return(list(chain_state = chain_state, logtarget = logtarget, logproposal = logproposal,
                accept_ratio = 1))
  }

  return(list(kernel = kernel))
}
