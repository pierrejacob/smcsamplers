### log-Gaussian Cox process
rm(list = ls())
library(MASS)
library(Rcpp)
library(RcppEigen)
library(temperingsmc)
library(spatstat)
library(tictoc)
library(mgcv)

### load in the pine forest data set
data(finpines)
### normalize the data to unit square
data_x <- (finpines$x + 5) / 10
data_y <- (finpines$y + 8) / 10
### specify the number of grids M
ngrid <- 10
grid <- seq(from = 0, to = 1, length.out = ngrid + 1)
dimension <- ngrid^2
### calculate the number of pine saplings in each grid cell
data_counts <- rep(0, dimension)
for (i in 1:ngrid){
  for (j in 1:ngrid){
    logical_y <- (data_x > grid[i]) * (data_x < grid[i+1])
    logical_x <- (data_y > grid[j]) * (data_y < grid[j+1])
    data_counts[(i-1)*ngrid + j] <- sum(logical_y * logical_x)
  }
}
### specify the parameters for the Gaussian process prior
parameter_sigmasq <- 1.91
parameter_beta <- 1/33
parameter_area <- 1 / dimension
parameter_mu <- log(126) - 0.5 * parameter_sigmasq
prior_mean <- rep(parameter_mu, dimension)
prior_cov <- matrix(nrow = dimension, ncol = dimension)
for(m in 1:dimension){
  for(n in 1:dimension){
    index_m <- c(floor((m - 1)/ngrid) + 1, ((m - 1)%%ngrid) + 1)
    index_n <- c(floor((n - 1)/ngrid) + 1, ((n - 1)%%ngrid) + 1)
    prior_cov[m, n] <- parameter_sigmasq * exp(-sqrt(sum((index_m - index_n)^2))/(ngrid * parameter_beta))
  }
}
prior_precision <- solve(prior_cov)
prior_precision_chol <- t(chol(prior_precision))
prior <- list()
prior$logdensity <- function(x) return(dmvnorm_cholesky_inverse(x, prior_mean, prior_precision_chol))
prior$rinit <- function(n) return(return(rmvn(n = n, mu = prior_mean, V = prior_cov)))
prior$gradlogdensity <- function(x){
  nparticles <- nrow(x)
  mu <- matrix(parameter_mu, nrow = nparticles, ncol = dimension)
  return(eigenMapMatMult(mu - x, prior_precision))
}

### compute metric tensor as in Girolami and Calderhead 2011 (Section 9)
### takes a few minutes to compute
metric <- list()
metric$tensor <- prior_precision
diag(metric$tensor) <- parameter_area * exp(prior_mean + 0.5 * diag(prior_cov)) + diag(prior_precision)
metric$inverse <- solve(metric$tensor)
cat("Inverse done")
metric$chol_inverse <- t(chol(metric$inverse))
cat("Chol done")
metric$inverse_chol_inverse <- solve(t(metric$chol_inverse))

### specify the target distribution
likelihood <- list()
likelihood$log <- function(x) return(coxprocess_loglikelihood(x, data_counts, parameter_area))
likelihood$gradlog <- function(x){
  nparticles <- nrow(x)
  return(matrix(rep(data_counts, nparticles), nrow = nparticles, byrow = T) - parameter_area * exp(x))
}
target_logdensity <- function(x) return(prior$logdensity(x) + likelihood$log(x))
target_gradlogdensity <- function(x) return(prior$gradlogdensity(x) + likelihood$gradlog(x))
target <- list(logdensity = target_logdensity, gradlogdensity = target_gradlogdensity)


### SMC algorithm
### set the proposal distribution to be the prior
proposal <- list(rinit = prior$rinit, logdensity = prior$logdensity, gradlogdensity = prior$gradlogdensity)
### we use the Remannian manifold HMC kernel for the MCMC rejuvenation step
mcmc_kernel <- list()
mcmc_kernel$choice <- "rm_hmc"
mcmc_kernel$parameters <- list(stepsize = 0.25, nsteps = 10, dimension = dimension, metric = metric)
### specify the number of MCMC moves for each intermediate distribution
mcmc_kernel$nmcmc <- 10
### specify the number of particles
nparticles <- 2^10
### specify the desingated effective sample size criterion
ess_criterion <- 0.5
tic("RM_HMC run time:")
smc_output <- run_smc(target, proposal, nparticles, mcmc_kernel = mcmc_kernel, ess_criterion)
toc()
### the log normalizing constant estimates
smc_output$log_ratio_normconst
