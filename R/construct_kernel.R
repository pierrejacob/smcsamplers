#' @rdname construct_kernel
#' @title Construct the MCMC rejuvenation kernel
#' @param target list with keys: \code{logdensity} evaluates log target density, \code{gradlogdensity} returns its gradient
#' @param proposal list with keys: \code{rinit} samples from proposal, \code{logdensity} evaluates log proposal density, \code{gradlogdensity} returns its gradient
#' @param mcmc_kernel list with keys: \code{choice} MCMC method, \code{parameters} associated tuning parameters, \code{nmcmc} number of MCMC moves for each intermediate distribution
#' @param gamma_current the current inverse temperature
#' @export
construct_kernel <- function(target, proposal, mcmc_kernel, gamma_current){
  ### Hamiltonian Monte Carlo
  if (mcmc_kernel$choice == "hmc"){
      return(hmc(target, proposal, mcmc_kernel$parameters, gamma_current)$kernel)
  }
  ### Gibbs sampler
  if (mcmc_kernel$choice == "gibbs"){
      return(gibbs(target, proposal, mcmc_kernel$parameters, gamma_current)$kernel)
  }
  ### Riemannian manifold Hamiltonian Monte Carlo
  if(mcmc_kernel$choice =="rm_hmc"){
    return(rm_hmc(target, proposal, mcmc_kernel$parameters, gamma_current)$kernel)
  }
  ### a sliced Gibbs sampler
  if(mcmc_kernel$choice =="slice_gibbs"){
    return(slice_gibbs(target, proposal, mcmc_kernel$conditional_distribution, gamma_current)$kernel)
  }
}
