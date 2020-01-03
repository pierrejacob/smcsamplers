### multivariate normal distribution
rm(list = ls())
library(temperingsmc)
library(mvtnorm)
library(mvnfast)
library(tictoc)


### specify the target distribution
dimension <- 50
rho <- 0.8
mean_target <- rep(0, dimension)
Sigma_target <- matrix(0, nrow = dimension, ncol = dimension)
for(i in 1:dimension){
  for(j in 1:dimension){
    Sigma_target[i, j] <- rho^abs(i - j)
  }
}
target <- get_mvnormal(dimension, mean_target, Sigma_target)


### SMC algorithm
### specify the proposal distribution
mean_proposal <- rep(0, dimension)
Sigma_proposal <- diag(sqrt(5), dimension, dimension)
proposal <- get_mvnormal(dimension, mean_proposal, Sigma_proposal)

### we use the HMC kernel for the MCMC rejuvenation step
### leapfrog generator: stepsize = 0.5, nsteps = 10
mcmc_kernel <- list()
mcmc_kernel$choice <- "hmc"
mcmc_kernel$parameters <- list(stepsize = 0.5, nsteps = 10)
### specify the number of MCMC moves for each intermediate distribution
mcmc_kernel$nmcmc <- 10
### specify the desingated effective sample size criterion
nparticles <- 2^10 + (dimension - 1) * 31
### specify the desingated effective sample size criterion
ess_criterion <- 0.8
tic("SMC runtime:")
smc_output <- run_smc(target, proposal, nparticles, mcmc_kernel, ess_criterion)
toc()
### the log normalizing constant estimates
### the truth is 0
smc_output$log_ratio_normconst


