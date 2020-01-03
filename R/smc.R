#' @rdname run_smc
#' @title Run sequential Monte Carlo sampler
#' @param target list with keys: \code{logdensity} evaluates log target density, \code{gradlogdensity} returns its gradient
#' @param proposal list with keys: \code{rinit} samples from proposal, \code{logdensity} evaluates log proposal density, \code{gradlogdensity} returns its gradient
#' @param nparticles number of particles
#' @param mcmc_kernel list with keys: \code{choice} MCMC method, \code{parameters} associated tuning parameters, \code{nmcmc} number of MCMC moves for each intermediate distribution
#' @param ess_criterion value between [0,1] which controls ESS adaptation
#' @export
run_smc <- function(target, proposal, nparticles, mcmc_kernel, ess_criterion){
  ### initialization
  ### sample from the initial distribution
  xparticles <- proposal$rinit(nparticles)
  xtrajectory <- list(xparticles)
  logtarget <- target$logdensity(xparticles)
  logproposal <- proposal$logdensity(xparticles)
  ### initialize the inverse temperature
  gamma_current <- 0
  temperature <- gamma_current
  ess <- 1
  ### initialize the log normalizing constant estimator
  log_ratio_normconst <- 0
  ### monitor the acceptance probability of the MCMC rejuvenation move
  accept_ratio <- numeric()
  while(gamma_current < 1){
    ### find the next inverse temperature
    logweights_IS <- logtarget - logproposal
    next_temp <- find_next_temp(gamma_current, ess_criterion, logweights_IS)
    ess <- c(ess, next_temp$ess)
    gamma_next <- next_temp$gamma_next
    temperature <- c(temperature, gamma_next)
    ### update the log normalizing constant estimator
    log_ratio_normconst <- log_ratio_normconst + next_temp$log_ratio_normconst
    normweights <- next_temp$normweights

    ### systematic resampling
    ancestors <- systematic_resampling(normweights, nparticles, runif(1))
    xparticles <- xparticles[ancestors, ]
    logtarget <- logtarget[ancestors]
    logproposal <- logproposal[ancestors]

    ### MCMC rejuvenation move
    transition_kernel <- construct_kernel(target, proposal, mcmc_kernel, gamma_next)
    for (i in 1:mcmc_kernel$nmcmc){
      if(mcmc_kernel$choice == "slice_gibbs"){
        rejuv_result <- transition_kernel(xparticles)
      }else{
        rejuv_result <- transition_kernel(xparticles, logtarget, logproposal)
      }
      xparticles <- rejuv_result$chain_state
      logtarget <- rejuv_result$logtarget
      logproposal <- rejuv_result$logproposal
      accept_ratio <- c(accept_ratio, rejuv_result$accept_ratio)
    }
    xtrajectory <- c(xtrajectory, list(xparticles))
    gamma_current <- gamma_next
  }

  return(list(xparticles = xparticles, log_ratio_normconst = log_ratio_normconst,
              xtrajectory = xtrajectory, temperature = temperature, ess = ess,
              accept_ratio = accept_ratio))
}

