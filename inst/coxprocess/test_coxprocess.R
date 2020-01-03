## This is an implementation of an adaptive SMC sampler
## on a tempered sequence bridging prior and posterior
## for the log-Gaussian Cox process example
## with slice sampling moves

rm(list = ls())
set.seed(1)

library(temperingsmc)
library(spatstat)
library(RcppTN)
### load in the pine forest data set
data(finpines)
### normalize data to unit square
data_x <- (finpines$x + 5) / 10
data_y <- (finpines$y + 8) / 10
# plot(data_x, data_y)

### specify a grid
ngrid <- 20
grid <- seq(from = 0, to = 1, length.out = ngrid + 1)
dimension <- ngrid^2
## might be convenient to introduce functions to map
## an index m in {1:ngrid}^2 to a pair (x,y) with x and y in 1:ngrid
m2xy <- function(m) c(floor((m - 1)/ngrid) + 1, ((m - 1)%%ngrid) + 1)
xy2m <- function(x,y) (x-1)*ngrid + y
# m2xy(10)
# m2xy(100)
# xy2m(m2xy(23)[1], m2xy(23)[2])
##

### calculate the number of pine saplings in each grid cell
## data_counts[(i-1)*ngrid + j] gives number of saplings in row i, column j
data_counts <- rep(0, dimension)
for (i in 1:ngrid){
  for (j in 1:ngrid){
    logical_y <- (data_x > grid[i]) * (data_x < grid[i+1])
    logical_x <- (data_y > grid[j]) * (data_y < grid[j+1])
    data_counts[xy2m(i,j)] <- sum(logical_y * logical_x)
  }
}

### specify the parameters for the Gaussian process prior
parameter_sigmasq <- 1.91
parameter_beta <- 1/33
parameter_area <- 1 / dimension
parameter_mu <- log(126) - 0.5 * parameter_sigmasq
prior_mean <- rep(parameter_mu, dimension)
## prior covariance: entry m,n equal to sigma^2 exp(-(1 / (M * beta)) * |m - n|)
prior_cov <- matrix(nrow = dimension, ncol = dimension)
for(m in 1:dimension){
  for(n in 1:dimension){
    xy_m <- m2xy(m)
    xy_n <- m2xy(n)
    prior_cov[m, n] <- parameter_sigmasq * exp(-sqrt(sum((xy_m - xy_n)^2))/(ngrid * parameter_beta))
  }
}
prior_cov_chol <- chol(prior_cov)
prior_precision <- solve(prior_cov)
prior_precision_chol <- t(chol(prior_precision))

## specify prior distribution
## also used as initial distribution
priordist <- list(logdensity = function(xs) temperingsmc:::dmvnorm_cholesky_inverse(xs, prior_mean, prior_precision_chol),
                    generate = function(n) temperingsmc:::rmvnorm_cholesky_(n, prior_mean, prior_cov_chol))

# prior <- list()
# prior$logdensity <- function(x) return(dmvnorm_cholesky_inverse(x, prior_mean, prior_precision_chol))
# prior$rinit <- function(n) return(return(rmvn(n = n, mu = prior_mean, V = prior_cov)))
# prior$gradlogdensity <- function(x){
#   nparticles <- nrow(x)
#   mu <- matrix(parameter_mu, nrow = nparticles, ncol = dimension)
#   return(eigenMapMatMult(mu - x, prior_precision))
# }

### specify the likelihood function
# likelihood <- list()
loglikelihood <- function(xs) coxprocess_loglikelihood(xs, data_counts, parameter_area)

# xs <- priordist$generate(2)
# loglikelihood(xs)

# likelihood$gradlog <- function(x) return(data_counts - parameter_area * exp(x))

# target_logdensity <- function(x) return(prior$logdensity(x) + likelihood$log(x))
# target_gradlogdensity <- function(x) return(prior$gradlogdensity(x) + likelihood$gradlog(x))
# target <- list(logdensity = target_logdensity, gradlogdensity = target_gradlogdensity)

## importance sampling from prior to posterior
# nparticles <- 1e6
# xparticles <- priordist$generate(nparticles)
# loglikes <- loglikelihood(xparticles)
# normalized_weights <- PET::normalize_weight(loglikes)
# ess <- 1/sum(normalized_weights$nw^2)
# normalized_weights$avew

### slice Gibbs sampling targeting the tempered posterior distribution
### pre-compute quantities for slice sampling
## slice_var[i] = sigma^2 - Sigma_{i,-i} Sigma_{-i,-i}^{-1} Sigma_{-i,i}
## slice_mean[,i] = Sigma_{i,-i} Sigma_{-i,-i}^{-1}
slice_mean <- matrix(0, nrow = dimension-1, ncol = dimension)
slice_var <- rep(0, dimension)
slice_sd <- rep(0, dimension)
for(d in 1:dimension){
  slice_mean[,d] <- eigenMapMatMult(matrix(prior_cov[d, -d], nrow = 1), solve(prior_cov[-d, -d]))
  slice_var[d] <- prior_cov[d, d] - eigenMapMatMult(matrix(slice_mean[,d], nrow = 1), prior_cov[-d, d])
  slice_sd[d] <- sqrt(slice_var[d])
}

### slice sampling move targeting tempered posterior (likelihood raised to power gamma)
## scans all the components of x-particles systematically
## and for each, sample a uniform given the x-particle and a new x-particle given the uniform
slice_move <- function(xparticles, gamma){
  nparticles <- dim(xparticles)[1]
  for (d in 1:dimension){
    # sample u given x
    us <- runif(nparticles, min = rep(0, nparticles), max = exp(-gamma * parameter_area * exp(xparticles[, d])))
    # sample x given u
    cond_sd <- rep(slice_sd[d], nparticles)
    cond_mean <- eigenMapMatMult(xparticles[, -d] - parameter_mu, slice_mean[,d,drop=F]) + parameter_mu
    # cond_mean <- eigenMapMatMult(xparticles[, -d] - matrix(parameter_mu, nrow = nparticles, ncol = dimension - 1),
    #                              matrix(slice_mean[d, ], ncol = 1)) + parameter_mu
    cond_mean <- cond_mean + gamma * data_counts[d] * slice_var[d]
    xparticles[,d] <- rtn(nparticles, .mean = cond_mean, .sd = cond_sd, .low = rep(-Inf, nparticles), .high = log(-1/(parameter_area * gamma) * log(us)))
  }
  return(xparticles)
}

# nsamples <- 2
# xparticles <- priordist$generate(nsamples)
# nmcmc <- 1000
# x1chain1 <- rep(0, nmcmc)
# x1chain2 <- rep(0, nmcmc)
# for (imcmc in 1:nmcmc){
#   xparticles <- slice_move(xparticles, 1)
#   x1chain1[imcmc] <- xparticles[1,1]
#   x1chain2[imcmc] <- xparticles[2,1]
# }
# matplot(cbind(x1chain1, x1chain2), type = "l")
# hist(x1chain1[(floor(nmcmc/10):nmcmc)], prob = TRUE)
# summary(x1chain1[(floor(nmcmc/10):nmcmc)]) # mean ~5.2


### adaptive SMC sampler
asmc_cox <- function(smctuning){
  ## start timer
  starttime <- Sys.time()
  attach(smctuning, warn.conflicts = FALSE)
  ## create particles: a list containing the locations of the particles "x"
  ## as well their target and initial log-density evaluations
  ## and 'n' for the number of particles, 'd' for the dimension
  particles <- list()
  particles$x <- priordist$generate(nparticles)
  particles$loglikelihood <- loglikelihood(particles$x)
  # particles$logprior <- priordist$logdensity(particles$x)
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
      logw <- (gamma - gamma_current) * particles$loglikelihood
      return(1/sum(PET::normalize_weight(logw)$nw^2))
    }
    if (ess_deltagamma(1) > ess_criterion * nparticles){
      gamma_next <- 1
    } else {
      search_results <- PET::search_gamma(gamma_current, ess_deltagamma, ess_criterion * nparticles)
      gamma_next <- search_results$x
    }
    ## now weight and resample
    logw <- (gamma_next - gamma_current) * particles$loglikelihood
    nw <- PET::normalize_weight(logw)
    log_ratio_normconst <- c(log_ratio_normconst, nw$avew)
    ## uses multinomial resampling, for variance estimation purposes
    ancestors <- PET::multinomial_resampling(nparticles, nw$nw)
    eves <- eves[ancestors]
    particles$x <- particles$x[ancestors,]
    particles$loglikelihood <- particles$loglikelihood[ancestors]
    # particles$logprior <- particles$logprior[ancestors]
    gamma_current <- gamma_next
    gammas <- c(gammas, gamma_current)
    ## MCMC move at current gamma
    for (imove in 1:nmoves){
      particles$x <- slice_move(particles$x, gamma_current)
    }
    particles$loglikelihood <- loglikelihood(particles$x)
  }
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(particles = particles, eves = eves, gammas = gammas, log_ratio_normconst = log_ratio_normconst, elapsedtime = elapsedtime))
}

smctuning <- list(nparticles = 1000, ess_criterion = 0.8, nmoves = 5)
asmcresults <- asmc_cox(smctuning)
asmcresults$gammas
sum(asmcresults$log_ratio_normconst)
asmcresults$elapsedtime

hist(asmcresults$particles$x[,1], prob = TRUE)
summary(asmcresults$particles$x[,1]) # mean ~5.2

# unique(asmcresults$eves)

