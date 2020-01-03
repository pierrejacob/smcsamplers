### WORK IN PROGRESS

## This is an implementation of an adaptive SMC sampler
## for the case of tempered distributions, with MALA moves
##
rm(list = ls())
set.seed(1)
##

## input:
## initial distribution pi_0 (how to sample from it, how to evaluate its density)
## target distribution pi (how to evaluate its density)

## let's define density functions such that
## - they can take as arguments all the particles in a N x d matrix where N is number of particles
## - they return values on log scale

dimension <- 12

## initial distribution and target distributions are specified as lists
## with 'logdensity' to compute log density at each particle
## for the initial distribution there's also a 'generate' function to start the SMC algorithm

## initial Normal(0, I)
initialdist_mean <- rep(0, dimension)
initialdist_variance <- diag(1, dimension, dimension)
initialdist_precision <- solve(initialdist_variance)
initialdist_cholvariance <- chol(initialdist_variance)
initialdist_cholinvvariance <- t(chol(solve(initialdist_variance)))
initialdist <- list(logdensity = function(xs) temperingsmc:::dmvnorm_cholesky_inverse(xs, initialdist_mean, initialdist_cholinvvariance),
                    generate = function(n) temperingsmc:::rmvnorm_cholesky_(n, initialdist_mean, initialdist_cholvariance))

## gradient of log density
initialdist$gradlogdensity <- function(xs) temperingsmc:::grad_dmvnorm_precision(xs, mean = initialdist_mean, initialdist_precision)
## test:
# initialdist$logdensity(initialdist$generate(2))

## target distribution Normal(mu, Sigma)
##  mu = c(0,1,0,1,....)
##  Sigma_ij = rho^{|i-j|}
rho <- 0.8
# targetdist_mean <- rep(1, dimension)
targetdist_mean <- rep(c(0,1), dimension/2)
targetdist_variance <- matrix(0, nrow = dimension, ncol = dimension)
for(i in 1:dimension) for(j in 1:dimension) targetdist_variance[i, j] <- rho^(abs(i - j))
targetdist_precision <- solve(targetdist_variance)
targetdist_cholvariance <- chol(targetdist_variance)
targetdist_cholinvvariance <- t(chol(solve(targetdist_variance)))
targetdist <- list(logdensity = function(xs) temperingsmc:::dmvnorm_cholesky_inverse(xs, targetdist_mean, targetdist_cholinvvariance))
targetdist$gradlogdensity <- function(xs) temperingsmc:::grad_dmvnorm_precision(xs, mean = targetdist_mean, targetdist_precision)

## MCMC
## the list 'mcmcmove' will contain functions to be used within the SMC sampler
## 'adapt' will take particles and produce tuning parameters
## 'sample' will perform one step of MCMC for each particle given as input, and
## leaving the target distribution pi_gamma invariant
mcmcmove <- list()
## take particles and produce tuning parameters
mcmcmove$adapt <- function(particles){
  tuning <- list(invmass = diag(diag(cov(particles$x))))
  tuning$mass <- solve(tuning$invmass)
  tuning$mass_chol <- chol(tuning$mass)
  tuning$mass_cholinv <- t(chol(solve(tuning$mass)))
  tuning$stepsize <- 0.1
  tuning$nsteps <- 1
  return(tuning)
}

mcmcmove$sample <- function(particles, targetdist, initialdist, gamma, tuning){
  ## generate momenta variables
  initial_momenta <- temperingsmc:::rmvnorm_cholesky_(particles$n, rep(0, particles$d), tuning$mass_chol)
  grad_ <- function(xs) (1-gamma) * initialdist$gradlogdensity(xs) + gamma * targetdist$gradlogdensity(xs)
  positions <- particles$x
  ##
  # leap frog integrator
  momenta <- initial_momenta + tuning$stepsize * grad_(positions) / 2
  for (step in 1:tuning$nsteps){
    positions <- positions + tuning$stepsize * momenta %*% tuning$invmass
    if (step != tuning$nsteps){
      momenta <- momenta + tuning$stepsize * grad_(positions)
    }
  }
  momenta <- momenta + tuning$stepsize * grad_(positions) / 2
  ## Now MH acceptance step
  xproposals_initial_pdfs <- initialdist$logdensity(positions)
  xproposals_target_pdfs <- targetdist$logdensity(positions)
  proposed_pdfs <- (1 - gamma) * xproposals_initial_pdfs + gamma * xproposals_target_pdfs
  current_pdfs  <- (1 - gamma) * particles$loginit + gamma * particles$logtarget
  mhratios <- proposed_pdfs - current_pdfs
  mhratios <- mhratios + temperingsmc:::dmvnorm_cholesky_inverse(momenta, rep(0, particles$d), tuning$mass_cholinv) -
    temperingsmc:::dmvnorm_cholesky_inverse(initial_momenta, rep(0, particles$d), tuning$mass_cholinv)
  if (any(is.na(mhratios))) mhratios[is.na(mhratios)] <- -Inf
  accepts <- log(runif(particles$n)) < mhratios
  particles$x[accepts,] <- positions[accepts,]
  particles$loginit[accepts] <- xproposals_initial_pdfs[accepts]
  particles$logtarget[accepts] <- xproposals_target_pdfs[accepts]
  return(particles)
}

# ### Try MCMC move
# mcmcmove$tuning <- list(invmass = diag(1, dimension, dimension))
# mcmcmove$tuning$mass <- solve(mcmcmove$tuning$invmass)
# mcmcmove$tuning$mass_chol <- chol(mcmcmove$tuning$mass)
# mcmcmove$tuning$mass_cholinv <- t(chol(solve(mcmcmove$tuning$mass)))
# mcmcmove$tuning$stepsize <- 0.1
# mcmcmove$tuning$nsteps <- 10
#
# chains <- list()
# chains$n <- 5
# chains$d <- dimension
# chains$x <- initialdist$generate(chains$n)
# chains$loginit <- initialdist$logdensity(chains$x)
# chains$logtarget <- targetdist$logdensity(chains$x)

# nmcmc <- 5000
# chains_history <- matrix(nrow = nmcmc, ncol = dimension)
# for (imcmc in 1:nmcmc){
#   chains <- mcmcmove$sample(particles = chains, targetdist, initialdist, gamma = 0.1, tuning = mcmcmove$tuning)
#   chains_history[imcmc,] <- chains$x[1,]
# }
# matplot(chains_history[1000:nmcmc,], type = "l")
#
# marginal_index <- 1
# hist(chains_history[1000:nmcmc,marginal_index], prob = TRUE, nclass = 50)
# curve(dnorm(x, mean = targetdist_mean[marginal_index], sd = sqrt(targetdist_variance[marginal_index,marginal_index])), add = T)
# cov(chains_history[1500:nmcmc,])

## adaptive SMC sampler
asmc <- function(targetdist, initialdist, smctuning, mcmcmove){
  starttime <- Sys.time()
  attach(smctuning, warn.conflicts = FALSE)
  particles <- list()
  particles$x <- initialdist$generate(nparticles)
  # xtrajectory <- list(xparticles)
  particles$logtarget <- targetdist$logdensity(particles$x)
  particles$loginit <- initialdist$logdensity(particles$x)
  particles$n <- nparticles
  particles$d <- dim(particles$x)[2]
  ### initialize the inverse temperature
  gamma_current <- 0
  gammas <- c(gamma_current)
  ### initialize the log normalizing constant estimator
  log_ratio_normconst <- c()
  ### monitor the acceptance probability of the MCMC rejuvenation move
  accept_ratio <- numeric()
  # while inv temperature is not one...
  while(gamma_current < 1){
    ### find the next inverse temperature
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
  }
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(particles = particles, gammas = gammas, log_ratio_normconst = log_ratio_normconst, elapsedtime = elapsedtime))
}

smctuning <- list()
smctuning$nparticles <- 2^10 + (dimension - 1) * 31

smctuning$ess_criterion <- 0.8
smctuning$nmoves <- dimension

asmcresults <- temperingsmc::asmc_tempering(targetdist, initialdist, smctuning, mcmcmove)
asmcresults$gammas
sum(asmcresults$log_ratio_normconst)
asmcresults$elapsedtime

marginal_index <- 2
hist(asmcresults$particles$x[,marginal_index], prob = TRUE, nclass = 50)
curve(dnorm(x, mean = targetdist_mean[marginal_index], sd = sqrt(targetdist_variance[marginal_index,marginal_index])), add = T)

## many repeats
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
nrep <- 40
asmc_ <- foreach(irep = 1:nrep) %dorng% asmc(targetdist, initialdist, smctuning, mcmcmove)
lognormcsts <- sapply(asmc_, function(x) sum(x$log_ratio_normconst))
hist(lognormcsts, prob = TRUE)
sd(lognormcsts)
