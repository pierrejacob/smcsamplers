# https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html#3_scaling_up_ode-based_models
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores ())
set.seed(3) # for reproductibility
stan_model_string <- "

functions {
  real[] sir(real t, real[] y, real[] theta,
             real[] x_r, int[] x_i) {

      real S = y[1];
      real I = y[2];
      real R = y[3];
      real N = x_i[1];

      real beta = theta[1];
      real gamma = theta[2];

      real dS_dt = -beta * I * S / N;
      real dI_dt =  beta * I * S / N - gamma * I;
      real dR_dt =  gamma * I;

      return {dS_dt, dI_dt, dR_dt};
  }
}

data {
  int<lower=1> n_days;
  real y0[3];
  real t0;
  real ts[n_days];
  int N;
  int cases[n_days];
}

transformed data {
  real x_r[0];
  int x_i[1] = { N };
}

parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> phi_inv;
}

transformed parameters{
  real y[n_days, 3];
  real phi = 1. / phi_inv;
  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
}

model {
  //priors
  beta ~ normal(2, 1);
  gamma ~ normal(0.4, 0.5);
  phi_inv ~ exponential(5);

  //sampling distribution
  //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people
  cases ~ neg_binomial_2(col(to_matrix(y), 2), phi);
}

generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  real pred_cases[n_days];
  pred_cases = neg_binomial_2_rng(col(to_matrix(y), 2), phi);
}

"

library(outbreaks)
library(tidyverse)
# time series of cases
cases <- influenza_england_1978_school$in_bed # Number of students in bed
write(stan_model_string, file = "sirmodel.stan")
sir_mod <- stan_model(file = "sirmodel.stan", verbose = F)

### total count
N <- 763
get_data_sir <- function(n_days, t0 = 0){
  # times
  t <- seq(0, n_days, by = 1)
  t0 = t0
  t <- t[-1]
  #initial conditions
  i0 <- 1
  s0 <- N - i0
  r0 <- 0
  y0 = c(S = s0, I = i0, R = r0)
  # data for Stan
  data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, cases = cases[1:n_days])
  return(data_sir)
}

library(smcsamplers)
graphsettings <- set_custom_theme()

## define functions for SMC samplers
get_sir_targetdist <- function(fit_){
  gradlogtarget_ <- function(x) grad_log_prob(fit_, x)
  targetdist <- function(xs){
    gradlogpdf <- matrix(NA, nrow = dim(xs)[1], ncol = dim(xs)[2])
    logpdf <- rep(NA, dim(xs)[2])
    for (ix in 1:ncol(xs)){
      grad_ <- try(gradlogtarget_(xs[,ix]), silent = T)
      if (inherits(grad_, "try-error")){
      } else {
        gradlogpdf[,ix] <- grad_
        logpdf[ix] <- attr(grad_, "log_prob")
      }
    }
    # apply(xs, 2, stan_gradlogtarget)
    return(list(gradlogpdf = gradlogpdf, logpdf = logpdf))
  }
  return(targetdist)
}

### start by going from prior to posterior given all data, via tempering of the likelihood
fit_sir_negbin <- sampling(sir_mod, data = get_data_sir(length(cases)), iter = 2, chains = 1)
targetdist <- get_sir_targetdist(fit_sir_negbin)

load("experiments/sirmodel/boardingschool-smc-partial.RData")
nsteps <- length(smc_partial_results[[1]])
xs <- (smc_partial_results[[1]])[[nsteps]]
postmean <- rowMeans(xs)
postcov <- cov(t(xs))

init <- get_mvnormal(postmean, postcov)

initdist_generate <- function(n) init$generate(n)
initdist <- init$eval

smctuning <- list(nparticles = 2^12, ess_criterion = 0.5, nmoves = 3)
smctuning$stepsize <- 0.1
smctuning$nleapfrog <- ceiling(1/smctuning$stepsize)
initparticles <- initdist_generate(smctuning$nparticles)
results <- asmc_hmc(smctuning, targetdist, initdist, initparticles)
nroots <- sapply(results$roots_history, function(x) length(unique(x)))
print(nroots)
results$ess_realized

t0s <- seq(from = -0.5, to = 0.8, length.out = 11)
smcresults_list <- list()
library(doParallel)
library(doRNG)
registerDoParallel(cores = 5)
nrep <- 5

for (it0 in seq_along(t0s)){
  t0 <- t0s[it0]
  cat("t0 =", t0, "\n")
  ## get compiled stan object
  fit_sir_negbin <- sampling(sir_mod, data = get_data_sir(length(cases), t0 = t0), iter = 2, chains = 1)
  ## get target distribution
  targetdist <- get_sir_targetdist(fit_sir_negbin)
  ## run SMC sampler from mv Normal to target
  smcresults <- foreach(irep = 1:nrep) %dorng% {
    initparticles <- initdist_generate(smctuning$nparticles)
    results <- asmc_hmc(smctuning, targetdist, initdist, initparticles)
  }
  smcresults_list[[it0]] <- smcresults
}
##

resultsdf <- data.frame()
for (it0 in seq_along(t0s)){
  t0 <- t0s[it0]
  sres <- smcresults_list[[it0]]
  zhats <- sapply(sres, function(x) sum(x$log_ratio_normconst))
  resultsdf <- rbind(resultsdf, data.frame(t0 = t0, zhat = zhats, irep = 1:nrep))
}

save(t0s, nrep, resultsdf, smctuning, file = "experiments/sirmodel/boardingschool-smc-marginal.RData")

head(resultsdf)
ggplot(resultsdf, aes(x = t0, y = zhat)) + geom_point() + xlab(expression(t[0])) + ylab("Z hat")



