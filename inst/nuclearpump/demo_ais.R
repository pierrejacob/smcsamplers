rm(list=ls())
library(ggplot2)

# prior specification
prior_gamma <- 0.01
prior_gamma_0 <- 4
prior_delta <- 1
prior_alpha <- 1.802
prior <- list()
prior$logdensity <- function(x) as.numeric(nuclear_logartificialprior(x))
prior$rinit <- function(n) nuclear_sample_artificialprior(n) 

# dataset
dataset_failures <- c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
dataset_times <- c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5)
nobservations <- length(dataset_failures) # number of pumps
dimension <- nobservations + 1 # dimension

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(nuclear_logartificiallikelihood(x)) 

# SMC settings
nparticles <- 2^8 # no. of particles
nsteps <- 10 # no. time steps
nmoves <- 5 # no. of MCMC moves
lambda <- seq(0, 1, length.out = nsteps+1) # inverse temperature schedule

# pre-allocate
logweights <- rep(0, nparticles)
ess <- rep(0, nsteps+1)
ess[1] <- nparticles
log_normconst <- rep(0, nsteps+1)

# initialize 
particles <- prior$rinit(nparticles)
betaa <- particles[,1]
theta <- particles[,2:dimension]

for (i in 1:nsteps){
  i
  # weight 
  current_loglikelihood <- likelihood$logdensity(particles)
  logweights <- logweights + (lambda[i+1] - lambda[i]) * current_loglikelihood
  maxlogweights <- max(logweights)
  weights <- exp(logweights - maxlogweights)
  normweights <- weights / sum(weights)
  
  # compute effective sample size
  ess[i+1] <- 1 / sum(normweights^2)
  
  # compute normalizing constant
  log_normconst[i+1] <- log(mean(weights)) + maxlogweights
  
  # MCMC
  for (m in 1:nmoves){
    # sample beta 
    betaa <- as.numeric(nuclear_sample_beta(rowSums(theta), lambda[i+1]))

    # sample theta
    theta <- nuclear_sample_theta(betaa, lambda[i+1])
    
  }
  particles <- cbind(betaa, theta)

}

# ess plot
ess.df <- data.frame(time = 0:nsteps, ess = ess * 100 / nparticles)
ggplot(ess.df, aes(x = time, y = ess)) + geom_line() + 
  labs(x = "time", y = "ESS%") + ylim(c(0, 100))

# normalizing constant plot
normconst.df <- data.frame(time = 0:nsteps, normconst = log_normconst)
ggplot() + geom_line(data = normconst.df, aes(x = time, y = normconst), colour = "blue") + 
  labs(x = "time", y = "log normalizing constant")
log_normconst[nsteps+1]
