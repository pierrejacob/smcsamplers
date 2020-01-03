## adaptive SMC sampler for tempered paths
## targetdist should be a list with 'logdensity', and possibly other functions required by 'mcmcmove$sample'
## initialdist should be a list with 'logdensity', 'generate', and possibly other functions required by 'mcmcmove$sample'
## smctuning should be a list with 'nparticles', 'nmoves', 'ess_criterion'
## mcmcmove should be a list with 'adapt', 'sample'
## Returns:
## * particles
## * eve indices
## * gammas
## * log_ratio_normconst
## * elapsedtime

#'@export
asmc_tempering <- function(targetdist, initialdist, smctuning, mcmcmove){
  ## start timer
  starttime <- Sys.time()
  attach(smctuning, warn.conflicts = FALSE)
  ## create particles: a list containing the locations of the particles "x"
  ## as well their target and initial log-density evaluations
  ## and 'n' for the number of particles, 'd' for the dimension
  particles <- list()
  particles$x <- initialdist$generate(nparticles)
  particles$logtarget <- targetdist$logdensity(particles$x)
  particles$loginit <- initialdist$logdensity(particles$x)
  particles$n <- nparticles
  particles$d <- dim(particles$x)[2]
  ## initialize the 'inverse temperature' at zero
  gamma_current <- 0
  ## store inverse temperatures
  gammas <- c(gamma_current)
  ## store log normalizing constant estimators at each step
  log_ratio_normconst <- c()
  ## eve indices
  eves <- 1:nparticles
  ## while inverse temperature is not one...
  while(gamma_current < 1){
    ### find next inverse temperature
    ess_deltagamma <- function(gamma){
      logw <- (gamma - gamma_current) * (particles$logtarget - particles$loginit)
      return(1/sum(PET::normalize_weight(logw)$nw^2))
    }
    if (ess_deltagamma(1) > ess_criterion * nparticles){
      gamma_next <- 1
    } else {
      search_results <- PET::search_gamma(gamma_current, ess_deltagamma, ess_criterion * nparticles)
      gamma_next <- search_results$x
    }
    ## now weight and resample
    logw <- (gamma_next - gamma_current) * (particles$logtarget - particles$loginit)
    nw <- PET::normalize_weight(logw)
    log_ratio_normconst <- c(log_ratio_normconst, nw$avew)
    ## uses multinomial resampling, for variance estimation purposes
    ancestors <- PET::multinomial_resampling(nparticles, nw$nw)
    eves <- eves[ancestors]
    particles$x <- particles$x[ancestors,]
    particles$logtarget <- particles$logtarget[ancestors]
    particles$loginit <- particles$loginit[ancestors]
    gamma_current <- gamma_next
    gammas <- c(gammas, gamma_current)
    ## MCMC move at current gamma
    mcmcmove$tuning <- mcmcmove$adapt(particles)
    for (imove in 1:nmoves){
      particles <- mcmcmove$sample(particles, targetdist, initialdist, gamma_current, mcmcmove$tuning)
    }
  }
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(particles = particles, eves = eves, gammas = gammas, log_ratio_normconst = log_ratio_normconst, elapsedtime = elapsedtime))
}


