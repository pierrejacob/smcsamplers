### log-Gaussian Cox process
rm(list = ls())
library(MASS)
library(Rcpp)
library(RcppEigen)
library(temperingsmc)
library(spatstat)
library(tictoc)
library(mgcv)
### the RcppTN package can be installed from source
### https://cran.r-project.org/web/packages/RcppTN/index.html
library(RcppTN)

### load in the pine forest data set
data(finpines)
### normalize data to unit square
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
### specify the target distribution
likelihood <- list()
likelihood$log <- function(x) return(coxprocess_loglikelihood(x, data_counts, parameter_area))
likelihood$gradlog <- function(x) return(data_counts - parameter_area * exp(x))
target_logdensity <- function(x) return(prior$logdensity(x) + likelihood$log(x))
target_gradlogdensity <- function(x) return(prior$gradlogdensity(x) + likelihood$gradlog(x))
target <- list(logdensity = target_logdensity, gradlogdensity = target_gradlogdensity)

### slice Gibbs sampling
### specify the parameters for slice sampling
slice_mean <- matrix(0, nrow = dimension, ncol = dimension - 1)
slice_var <- rep(0, dimension)
for(i in 1:dimension){
  slice_mean[i, ] <- eigenMapMatMult(matrix(prior_cov[i, -i], nrow = 1), solve(prior_cov[-i, -i]))
  slice_var[i] <- prior_cov[i, i] - eigenMapMatMult(matrix(slice_mean[i, ], nrow = 1), prior_cov[-i, i])
}
### conditional distribution of x given u
### truncated normal distribution
cond_x <- function(x, u, d, gamma){
  nparticles <- nrow(x)
  cond_var <- rep(slice_var[d], nparticles)
  cond_mean <- eigenMapMatMult(x[, -d] - matrix(parameter_mu, nrow = nparticles, ncol = dimension - 1),
                               matrix(slice_mean[d, ], ncol = 1)) + parameter_mu
  cond_mean <- cond_mean + gamma * data_counts[d] * slice_var[d]
  return(rtn(nparticles, .mean = cond_mean, .sd = sqrt(cond_var), .low = rep(-Inf, nparticles), .high = log(-1/(parameter_area * gamma) * log(u))))
}
### conditional distribution of u given x
### uniform distribution
cond_u <- function(x, d, gamma){
  nparticles <- nrow(x)
  return(runif(nparticles, min = rep(0, nparticles), max = exp(-gamma * parameter_area * exp(x[, d]))))
}

### SMC algorithm
### set the proposal distribution to be the prior
proposal <- prior
### we use the slice Gibbs kernel for the MCMC rejuvenation step
mcmc_kernel <- list()
mcmc_kernel$choice <- "slice_gibbs"
mcmc_kernel$conditional_distribution <- list(conditional_x = cond_x, conditional_u = cond_u)
### specify the number of MCMC moves for each intermediate distribution
mcmc_kernel$nmcmc <- 10
### specify the number of particles
nparticles <- 2^10
### specify the desingated effective sample size criterion
ess_criterion <- 0.5
tic("slice gibbs run time:")
smc_output <- run_smc(target, proposal, nparticles, mcmc_kernel = mcmc_kernel, ess_criterion)
toc()
### the log normalizing constant estimates
smc_output$log_ratio_normconst

