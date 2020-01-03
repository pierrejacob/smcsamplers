#'@rdname rm_hmc
#'@title Run Riemann manifold Hamiltonian Monte Carlo
#'@param target list with keys: \code{logdensity} evaluates log target density, \code{gradlogdensity} returns its gradient
#'@param proposal list with keys: \code{rinit} samples from proposal, \code{logdensity} evaluates log proposal density, \code{gradlogdensity} returns its gradient
#'@param rm_hmc_parameters list with keys: \code{stepsize} step size in the leap-frog integrator, \code{nsteps} number of leapfrog steps
#'@param gamma_current the current inverse temperature
#'@export
rm_hmc <- function(target, proposal, rm_hmc_parameters, gamma_current){
  stepsize <- rm_hmc_parameters$stepsize
  nsteps <- rm_hmc_parameters$nsteps
  dimension <- rm_hmc_parameters$dimension
  metric <- rm_hmc_parameters$metric
  gradient <- function(x) gamma_current * target$gradlogdensity(x) + (1 - gamma_current) * proposal$gradlogdensity(x)

  ### leap frog integrator
  leapfrog <- function(x, v){
    v <- v + stepsize * gradient(x) / 2
    for (step in 1:nsteps){
      x <- x + stepsize * (v%*%t(metric$inverse))
      if (step != nsteps){
        v <- v + stepsize * gradient(x)
      }
    }
    v <- v + stepsize * gradient(x) / 2
    ### we could negate the momentum but we don't use it here
    return(list(x = x, v = v))
  }

  ### construct the RM-HMC kernel
  kernel <- function(chain_state, logtarget_current, logproposal_current){
    ### chain_state: N x d matrix; logtarget_current: N x 1 vector; logproposal_current: N x 1 vector
    nparticles <- nrow(chain_state)
    ### refresh the velocity
    current_v <- matrix(rnorm(dimension * nparticles), ncol = dimension)%*%metric$inverse_chol_inverse
    ### simulate the Hamiltonian dynamic
    leapfrog_result <- leapfrog(chain_state, current_v)
    proposed_v <- - leapfrog_result$v
    proposed_x <- leapfrog_result$x

    ### the Metropolis acceptance and rejection step
    logtarget_proposed <- target$logdensity(proposed_x)
    logproposal_proposed <- proposal$logdensity(proposed_x)
    loglik_proposed <- gamma_current * logtarget_proposed + (1 - gamma_current) * logproposal_proposed
    loglik_current <- gamma_current * logtarget_current + (1 - gamma_current) * logproposal_current
    accept_ratio <- logtarget_proposed - logtarget_current
    ### the acceptance ratio also features the "kinetic energy" term of the extended target
    accept_ratio <- accept_ratio + 0.5 * rowSums((current_v %*% metric$inverse) * current_v) -
                                   0.5 * rowSums((current_v %*% metric$inverse) * current_v)
    accept_logical <- log(runif(nparticles)) < accept_ratio
    chain_state[accept_logical, ] <- proposed_x[accept_logical, ]
    logtarget_current[accept_logical] <- logtarget_proposed[accept_logical]
    logproposal_current[accept_logical] <- logproposal_proposed[accept_logical]

    return(list(chain_state = chain_state, logtarget = logtarget_current, logproposal = logproposal_current,
                accept_ratio = mean(pmin(1, exp(accept_ratio)))))
  }

  return(list(kernel = kernel))
}
