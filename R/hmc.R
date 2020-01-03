#' @rdname hmc
#' @title Run Hamiltonian Monte Carlo
#' @param target list with keys: \code{logdensity} evaluates log target density, \code{gradlogdensity} returns its gradient
#' @param proposal list with keys: \code{rinit} samples from proposal, \code{logdensity} evaluates log proposal density, \code{gradlogdensity} returns its gradient
#' @param hmc_parameters list with keys: \code{stepsize} leapfrog step size, \code{nsteps}, number of leapfrog steps
#' @param gamma_current the current inverse temperature
#' @export
hmc <- function(target, proposal, hmc_parameters, gamma_current){
  ### leapfrog step size
  stepsize <- hmc_parameters$stepsize
  ### number of leapfrog steps
  nsteps <- hmc_parameters$nsteps
  ### the gradient of the log density of the current intermediate distribution
  gradient <- function(x) gamma_current * target$gradlogdensity(x) + (1 - gamma_current) * proposal$gradlogdensity(x)

  ### leap frog integrator
  leapfrog <- function(x, v){
    ### x: N x d matrix; v: N x d matrix
    ### N: the number of particles
    ### d: the dimension of the problem
    v <- v + stepsize * gradient(x) / 2
    for (step in 1:nsteps){
      x <- x + stepsize * v
      if (step != nsteps){
        v <- v + stepsize * gradient(x)
      }
    }
    v <- v + stepsize * gradient(x) / 2
    ### we could negate the momentum but we don't use it here
    return(list(x = x, v = v))
  }

  ### construct the Hamiltonian Monte Carlo kernel
  kernel <- function(chain_state, logtarget_current, logproposal_current) {
    ### chain_state: N x d matrix; logtarget_current: N x 1 vector; logproposal_current: N x 1 vector
    nparticles <- nrow(chain_state)
    dimension <- ncol(chain_state)
    ### refresh the velocity
    current_v <- matrix(rnorm(dimension * nparticles), ncol = dimension)
    ### simulate the Hamiltonian dynamic
    leapfrog_result <- leapfrog(chain_state, current_v)
    proposed_v <- - leapfrog_result$v
    proposed_x <- leapfrog_result$x

    ### the Metropolis acceptance and rejection step
    logtarget_proposed <- target$logdensity(proposed_x)
    logproposal_proposed <- proposal$logdensity(proposed_x)
    loglik_proposed <- gamma_current * logtarget_proposed + (1 - gamma_current) * logproposal_proposed
    loglik_current <- gamma_current * logtarget_current + (1 - gamma_current) * logproposal_current
    ### calculate the log acceptance probability
    accept_ratio <- loglik_proposed - loglik_current + rowSums(current_v^2) / 2 - rowSums(proposed_v^2) / 2
    accept_logical <- (log(runif(nparticles)) < accept_ratio)
    accept_logical[is.na(accept_logical)] <- FALSE
    chain_state[accept_logical, ] <- proposed_x[accept_logical, ]
    logtarget_current[accept_logical] <- logtarget_proposed[accept_logical]
    logproposal_current[accept_logical] <- logproposal_proposed[accept_logical]

    return(list(chain_state = chain_state, logtarget = logtarget_current, logproposal = logproposal_current,
                accept_ratio = mean(pmin(1, exp(accept_ratio)))))
  }
  return(list(kernel = kernel))
}
