## This is an implementation of an adaptive SMC sampler
## for the case of tempered distributions, with RWMH moves
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
initialdist_cholvariance <- chol(initialdist_variance)
initialdist_cholinvvariance <- t(chol(solve(initialdist_variance)))
initialdist <- list(logdensity = function(xs) temperingsmc:::dmvnorm_cholesky_inverse(xs, initialdist_mean, initialdist_cholinvvariance),
                    generate = function(n) temperingsmc:::rmvnorm_cholesky_(n, initialdist_mean, initialdist_cholvariance))

## test:
# initialdist$logdensity(initialdist$generate(2))

## target distribution Normal(mu, Sigma)
##  mu = c(0,1,0,1,....)
##  Sigma_ij = rho^{|i-j|}
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

# ### Try MCMC move
# chains <- list()
# chains$n <- 2
# chains$d <- dimension
# chains$x <- initialdist$generate(chains$n)
# chains$loginit <- initialdist$logdensity(chains$x)
# chains$logtarget <- targetdist$logdensity(chains$x)
# mcmcmove$tuning <- list(chol_proposal_covariance = targetdist_cholvariance)
# nmcmc <- 100000
# chains_history <- matrix(nrow = nmcmc, ncol = dimension)
# for (imcmc in 1:nmcmc){
#   chains <- mcmcmove$sample(chains, targetdist, initialdist, 1, tuning = mcmcmove$tuning)
#   chains_history[imcmc,] <- chains$x[2,]
# }
# matplot(chains_history, type = "l")
# marginal_index <- 2
# hist(chains_history[500:nmcmc,marginal_index], prob = TRUE, nclass = 50)
# curve(dnorm(x, mean = targetdist_mean[marginal_index], sd = sqrt(targetdist_variance[marginal_index,marginal_index])), add = T)
# cov(chains_history[500:nmcmc,])


smctuning <- list()
smctuning$nparticles <- 2^10 + (dimension - 1) * 31

smctuning$ess_criterion <- 0.6
smctuning$nmoves <- dimension

asmcresults <- asmc_tempering(targetdist, initialdist, smctuning, mcmcmove)
asmcresults$gammas
sum(asmcresults$log_ratio_normconst)
asmcresults$elapsedtime

marginal_index <- 2
hist(asmcresults$particles$x[,marginal_index], prob = TRUE, nclass = 50, main = "marginal approximation", xlab = "x")
curve(dnorm(x, mean = targetdist_mean[marginal_index], sd = sqrt(targetdist_variance[marginal_index,marginal_index])), add = T)

## many repeats
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
nrep <- 40
asmc_ <- foreach(irep = 1:nrep) %dorng% asmc_tempering(targetdist, initialdist, smctuning, mcmcmove)
lognormcsts <- sapply(asmc_, function(x) sum(x$log_ratio_normconst))
hist(lognormcsts, prob = TRUE, main = "log norm constants")
sd(lognormcsts)
