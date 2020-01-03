set.seed(1)
dimension <- 60

## initial distribution and target distributions are specified as lists
## with 'logdensity' to compute log density at each particle
## for the initial distribution there's also a 'generate' function to start the SMC algorithm

## initial Normal(0, I)
initialdist_mean <- rep(0, dimension)
initialdist_variance <- diag(1, dimension, dimension)
initialdist_cholvariance <- chol(initialdist_variance)
initialdist_cholinvvariance <- t(chol(solve(initialdist_variance)))
initialdist <- list(logdensity = function(xs) temperingsmc:::dmvnorm_cholesky_inverse(xs, initialdist_mean, initialdist_cholinvvariance),
                    generate = function(n) temperingsmc:::rmvnorm_cholesky_(n, initialdist_mean, initialdist_cholvariance))
rho <- 0.8
targetdist_mean <- rep(c(0,1), dimension/2)
targetdist_variance <- matrix(0, nrow = dimension, ncol = dimension)
for(i in 1:dimension) for(j in 1:dimension) targetdist_variance[i, j] <- rho^(abs(i - j))
targetdist_cholvariance <- chol(targetdist_variance)
targetdist_cholinvvariance <- t(chol(solve(targetdist_variance)))
targetdist <- list(logdensity = function(xs) temperingsmc:::dmvnorm_cholesky_inverse(xs, targetdist_mean, targetdist_cholinvvariance))

## MCMC
## the list 'mcmcmove' will contain functions to be used within the SMC sampler
## 'adapt' will take particles and produce tuning parameters
## 'sample' will perform one step of MCMC for each particle given as input, and
## leaving the target distribution pi_gamma invariant
mcmcmove <- list()
## take particles and produce tuning parameters
mcmcmove$adapt <- function(particles){
  tuning <- list(proposal_covariance = cov(particles$x) / particles$d)
  tuning$chol_proposal_covariance <- chol(tuning$proposal_covariance)
  return(tuning)
}

## random walk MH kernel
## takes particles, target, initial distributions, gamma, and tuning parameters
## and output new particles
mcmcmove$sample <- function(particles, targetdist, initialdist, gamma, tuning){
  ## generate proposals
  xproposals <- particles$x + temperingsmc:::rmvnorm_cholesky_(particles$n, rep(0, particles$d), tuning$chol_proposal_covariance)
  ## evaluate pdf of target and initial at proposed points
  xproposals_initial_pdfs <- initialdist$logdensity(xproposals)
  xproposals_target_pdfs <- targetdist$logdensity(xproposals)
  ## pdf of current target
  xproposals_current_target_pdfs <- (1 - gamma) * xproposals_initial_pdfs + gamma * xproposals_target_pdfs
  current_target_pdfs <- (1 - gamma) * particles$loginit + gamma * particles$logtarget
  ## MH acceptance ratio on log scale
  mhratio <- (xproposals_current_target_pdfs - current_target_pdfs)
  ## decisions to accept proposed points
  accepts <- log(runif(particles$n)) < mhratio
  ## replacement of current particles by accepted proposed points
  particles$x[accepts,] <- xproposals[accepts,]
  particles$loginit[accepts] <- xproposals_initial_pdfs[accepts]
  particles$logtarget[accepts] <- xproposals_target_pdfs[accepts]
  return(particles)
}


##
library(rhdf5)
hdf5path <- "~/"
hdf5name <- "test_hdf5"
hdf5name <- paste(hdf5path, hdf5name, ".h5", sep="")
file.remove(hdf5name)

h5createFile(hdf5name)
h5createGroup(hdf5name, "history")
h5ls(hdf5name)

for (i in 1:10){
  samples <- temperingsmc:::rmvnorm_(nsamples, mean_, cov_)
  h5write(samples, file = hdf5name, name=paste0("history/samples", i))
}

#
h5ls(hdf5name)

h5ls(hdf5name)
h5read(hdf5name, "history/samples1")
h5read(hdf5name, "history/samples1", index = list(1:3,NULL))
h5read(hdf5name, "history/samples1", index = list(1:500,1:2))

h5closeAll()
#
# library(bench) #uses utils::Rprofmem
# gc()
# mark(x <- temperingsmc:::rmvnorm_(1e5, mean_, cov_))[,"mem_alloc"]
# #84MB

## now SMC storing everything
asmc_tempering_storeall <- function(targetdist, initialdist, smctuning, mcmcmove){
  ## start timer
  starttime <- Sys.time()
  attach(smctuning, warn.conflicts = FALSE)
  ##
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
  nsteps <- 1
  particles_history <- list()
  particles_history[[nsteps]] <- particles$x
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
    ancestors <- PET::ssp_resampling(nparticles, nw$nw)
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
    nsteps <- nsteps + 1
    particles_history[[nsteps]] <- particles$x
  }
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(particles_history = particles_history, particles = particles, gammas = gammas, log_ratio_normconst = log_ratio_normconst, elapsedtime = elapsedtime))
}

asmc_tempering_hdf5 <- function(targetdist, initialdist, smctuning, mcmcmove, hdf5file){
  ## start timer
  starttime <- Sys.time()
  attach(smctuning, warn.conflicts = FALSE)
  ##
  h5createFile(hdf5file)
  h5createGroup(hdf5file, "history")
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
  nsteps <- 1
  h5write(particles$x, file = hdf5file, name=paste0("history/particles", nsteps))
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
    ancestors <- PET::ssp_resampling(nparticles, nw$nw)
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
    nsteps <- nsteps + 1
    h5write(particles$x, file = hdf5file, name=paste0("history/particles_step", nsteps))
  }
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  h5closeAll()
  return(list(particles = particles, gammas = gammas, log_ratio_normconst = log_ratio_normconst, elapsedtime = elapsedtime))
}


smctuning <- list()
smctuning$nparticles <- 2^10 + (dimension - 1) * 31
smctuning$ess_criterion <- 0.99
smctuning$nmoves <- dimension

# file = "~/test.h5"
# file.remove(file)
# asmcresults <- asmc_tempering_hdf5(targetdist, initialdist, smctuning, mcmcmove, file)
# asmcresults$gammas
# asmcresults2 <- asmc_tempering_storeall(targetdist, initialdist, smctuning, mcmcmove)

library(doParallel)
library(doRNG)
registerDoParallel(cores = 10)
testxx <- foreach(irep = 1:10) %dorng% {
  file <- paste0("~/test", irep, ".h5")
  file.remove(file)
  asmcresults <- asmc_tempering_hdf5(targetdist, initialdist, smctuning, mcmcmove, file)
}

testyy <- foreach(irep = 1:10) %dorng% {
  asmcresults <- asmc_tempering_storeall(targetdist, initialdist, smctuning, mcmcmove)
}

# library(bench) #uses utils::Rprofmem
# file.remove(file)
# gc()
# mark(asmcresults <- asmc_tempering_hdf5(targetdist, initialdist, smctuning, mcmcmove, file))[,"mem_alloc"]
# h5ls(file)
# gc()
# file.remove(file)
# Rprof("Rprof.out", memory.profiling=TRUE, interval = runif(1,.005,0.02))
# asmcresults <- asmc_tempering_hdf5(targetdist, initialdist, smctuning, mcmcmove, file)
# Rprof(NULL)
# max(summaryRprof("Rprof.out", memory="both")$by.total$mem.total)
