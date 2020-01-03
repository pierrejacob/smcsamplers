## This is an implementation of an adaptive SMC sampler
## for the case of tempered distributions, with HMC moves
##
rm(list = ls())
library(doParallel)
library(doRNG)
library(temperingsmc)
registerDoParallel(cores = detectCores()-2)

set.seed(1)
##

## input:
## initial distribution pi_0 (how to sample from it, how to evaluate its density)
## target distribution pi (how to evaluate its density)

## let's define density functions such that
## - they can take as arguments all the particles in a N x d matrix where N is number of particles
## - they return values on log scale

dimension <- 4

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
  tuning$nsteps <- 10
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

smctuning <- list()
smctuning$nparticles <- 2^10 + (dimension - 1) * 31

smctuning$ess_criterion <- 0.8
smctuning$nmoves <- dimension

asmcresults <- temperingsmc::asmc_tempering(targetdist, initialdist, smctuning, mcmcmove)
gammas <- asmcresults$gammas
# sum(asmcresults$log_ratio_normconst)
# asmcresults$elapsedtime
# marginal_index <- 2
# hist(asmcresults$particles$x[,marginal_index], prob = TRUE, nclass = 50)
# curve(dnorm(x, mean = targetdist_mean[marginal_index], sd = sqrt(targetdist_variance[marginal_index,marginal_index])), add = T)

##
## non-adaptive SMC
smc_tempering <- function(targetdist, initialdist, smctuning, mcmcmove){
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
  igamma <- 1
  ## store log normalizing constant estimators at each step
  log_ratio_normconst <- c()
  ## eve indices
  eves <- 1:particles$n
  ## while inverse temperature is not one...
  for (igamma in 2:length(gammas)){
    # weight
    gamma_next <- gammas[igamma]
    logw <- (gamma_next - gamma_current) * (particles$logtarget - particles$loginit)
    nw <- PET::normalize_weight(logw)
    log_ratio_normconst <- c(log_ratio_normconst, nw$avew)
    # resample
    # ancestors <- PET::ssp_resampling(nparticles, nw$nw)
    ancestors <- PET::multinomial_resampling(nparticles, nw$nw)
    particles$x <- particles$x[ancestors,]
    particles$logtarget <- particles$logtarget[ancestors]
    particles$loginit <- particles$loginit[ancestors]
    eves <- eves[ancestors]
    gamma_current <- gamma_next
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

smctuning$nparticles <- 2^12
smctuning$gammas <- asmcresults$gammas
smcresults <- smc_tempering(targetdist, initialdist, smctuning, mcmcmove)

length(unique(smcresults$eves))

np1 <- length(smcresults$gammas)
## compute variance estimator as in lee & whiteley
## test function: phi maps any x to 1
prod_ <- (smctuning$nparticles/(smctuning$nparticles-1))^(np1)
sum_ <- 0
for (i in 1:smctuning$nparticles){
  sameeve <- which(smcresults$eves == i)
  sum_ <- sum_ + length(sameeve)^2
}
sum_/(smctuning$nparticles^2)

sums_ <- rep(0, smctuning$nparticles)
for (i in 1:smctuning$nparticles){
  sums_[smcresults$eves[i]] <- sums_[smcresults$eves[i]] + 1/smctuning$nparticles
}
sum(sums_^2)

var_estim <- 1 * (1 - prod_) + prod_ * sum_/(smctuning$nparticles^2)
var_estim
sqrt(var_estim)

compute_var_estim <- function(smcresults){
  nparticles <- length(smcresults$eves)
  np1 <- length(smcresults$gammas)
  prod_ <- (nparticles/(nparticles-1))^(np1)
  sums_ <- rep(0, nparticles)
  for (i in 1:nparticles){
    sums_[smcresults$eves[i]] <- sums_[smcresults$eves[i]] + 1/nparticles
  }
  var_estim <- 1 * (1 - prod_) + prod_ * sum(sums_^2)
  return(var_estim)
}
compute_var_estim(smcresults)

# Ns <- c(1e3, 2e3, 3e3, 4e3)
# nrep <- 200
# var_indeprepeats <- c()
# aggreg_var_estim <- c()
# singlerun_var_estim <- c()
# for (N in Ns){
#   print(N)
#   smctuning$nparticles <- N
#   smc_ <- foreach(irep = 1:nrep) %dorng% smc_tempering(targetdist, initialdist, smctuning, mcmcmove)
#   lognormcsts <- sapply(smc_, function(x) sum(x$log_ratio_normconst))
#   var_indeprepeats <- c(var_indeprepeats, var(lognormcsts))
#   var_estims <- sapply(smc_, function(x) compute_var_estim(x))
#   M <- mean(exp(lognormcsts))
#   aggreg_var_estim <- c(aggreg_var_estim, mean((exp(lognormcsts)/M)^2 * var_estims))
#   singlerun_var_estim <- c(singlerun_var_estim, var_estims[1])
# }
#
# var_indeprepeats
# aggreg_var_estim
# singlerun_var_estim
#
# plot(Ns, var_indeprepeats, ylim = c(0, 2.5e-3))
# lines(Ns, singlerun_var_estim)
# lines(Ns, aggreg_var_estim, col = 'red', lty = 2)

## what about adaptive SMC samplers?
smctuning$nparticles <- 2^13
asmcresults <- temperingsmc::asmc_tempering(targetdist, initialdist, smctuning, mcmcmove)
sqrt(compute_var_estim(asmcresults))
nrep <- 100
asmc_ <- foreach(irep = 1:nrep) %dorng% temperingsmc::asmc_tempering(targetdist, initialdist, smctuning, mcmcmove)
lognormcsts <- sapply(asmc_, function(x) sum(x$log_ratio_normconst))
sd_estims <- sapply(asmc_, function(x) sqrt(compute_var_estim(x)))
sd(lognormcsts)
summary(sd_estims)
## seems to work too!

