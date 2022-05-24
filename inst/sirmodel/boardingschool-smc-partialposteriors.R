# https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html#3_scaling_up_ode-based_models
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)
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
# time series of cases
cases <- influenza_england_1978_school$in_bed # Number of students in bed
write(stan_model_string, file = "sirmodel.stan")
sir_mod <- stan_model(file = "sirmodel.stan", verbose = F)


### total count
N <- 763
get_data_sir <- function(n_days){
  # times
  t <- seq(0, n_days, by = 1)
  t0 = 0
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

## start by doing MCMC targeting posterior distribution given 3 observations
fit_ <- sampling(sir_mod, data = get_data_sir(3), iter = 2, chains = 1)
library(smcsamplers)
graphsettings <- set_custom_theme()

## define target distribution
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
##
targetdist <- get_sir_targetdist(fit_)
# targetdist(matrix(c(-1,1,-1), ncol = 1))
##
rwmh_ <- function(init, Sigma, niterations){
  current <- init
  current_target <- targetdist(current)
  history <- matrix(ncol = niterations, nrow = 3)
  history[,1] <- current
  ##
  for (iteration in 2:niterations){
    if (iteration %% 100 == 1){
      cat("iteration:", iteration, "/", niterations, "\n")
    }
    proposal <- matrix(smcsamplers:::rmvnorm_(1, mean = t(current), Sigma)[1,], ncol = 1)
    proposal_target <- targetdist(proposal)
    u_ <- runif(1)
    mhratio <- (proposal_target$logpdf - current_target$logpdf)
    if (is.na(mhratio)) mhratio <- -Inf
    if (log(u_) < mhratio){
      current <- proposal
      current_target <- proposal_target
    }
    history[,iteration] <- current
  }
  return(history)
}
init <- matrix(c(-1,1,-1), ncol = 1)
Sigma <- diag(c(0.01, 0.01, 0.1))
niterations <- 2000
chain <- rwmh_(init, Sigma, niterations)
Sigma <- cov(t(chain))
chain <- rwmh_(chain[,niterations,drop=F], Sigma, niterations)
Sigma <- cov(t(chain))
chain <- rwmh_(chain[,niterations,drop=F], Sigma, niterations)
cat('acceptance rate:', 100*(1-mean(abs(diff(chain[1,]))<1e-10)), "%\n")
matplot(t(chain), type = 'l')

# ## we can compare the chain with what Stan gives
# load("experiments/sirmodel/boardingschool-stan-ndays3.RData")
# par(mfrow = c(3,1))
# hist(chain[1,], prob = T, nclass = 50)
# hist(log(standf$gamma), add = T, prob = T, col = rgb(1,0,0,0.5), nclass = 50)
# hist(chain[2,], prob = T, nclass = 50)
# hist(log(standf$beta), add = T, prob = T, col = rgb(1,0,0,0.5), nclass = 50)
# hist(chain[3,], prob = T, nclass = 50)
# hist(log(standf$phi_inv), add = T, prob = T, col = rgb(1,0,0,0.5), nclass = 50)
## close agreement, it seems

## now SMC sampler approximation
initdist_ <- get_mvnormal(rowMeans(chain), cov(t(chain)))

smctuning <- list(nparticles = 2^8, ess_criterion = 0.5, nmoves = 2)
smctuning$stepsize <- 0.1
smctuning$nleapfrog <- ceiling(1/smctuning$stepsize)
initparticles <- initdist_$generate(smctuning$nparticles)
results <- asmc_hmc(smctuning, targetdist, initdist_$eval, initparticles)
nroots <- sapply(results$roots_history, function(x) length(unique(x)))
print(nroots)
## get information about how the MCMC steps performed
info_mcmc_df <- lapply(2:length(results$lambdas), function(ilambda) {
  info_mcmc <- results$infos_mcmc[[ilambda-1]]
  df_ <- data.frame(ilambda = ilambda, imove = 1:smctuning$nmoves,
                    ar = sapply(info_mcmc, function(x) x$ar),
                    sqjd = sapply(info_mcmc, function(x) x$sqjd))
}) %>% bind_rows()
info_mcmc_df
##

## it seems to work so
## we will follow a path of distributions
## from the initial Gaussian approximation to the partial posteriors from 3 observations onwards

initmean <- c(-0.97, 0.49, -2.6)
initcov <- structure(c(1.59, 0.25, 0.16,
            0.25, 0.06, 0.04, 0.16,
            0.04, 1.08), .Dim = c(3L, 3L))
initdist_ <- get_mvnormal(initmean, initcov)
smctuning <- list(nparticles = 2^9, ess_criterion = 0.5, nmoves = 2)
smctuning$stepsize <- 0.1
smctuning$nleapfrog <- ceiling(1/smctuning$stepsize)
nrep <- 5
registerDoParallel(cores = 6)

smc_partial_results <- foreach(irep = 1:nrep) %dorng% {
  initparticles <- initdist_$generate(smctuning$nparticles)
  fit_ <- sampling(sir_mod, data = get_data_sir(3), iter = 2, chains = 1)
  targetdist <- get_sir_targetdist(fit_)
  results <- asmc_hmc(smctuning, targetdist, initdist_$eval, initparticles)
  xparticles <- results$particles$x
  xhistory <- list()
  xhistory[[1]] <- xparticles
  for (nobs in 4:length(cases)){
    cat("observation", nobs, "\n")
    ### start by going from prior to posterior given all data, via tempering of the likelihood
    # fit_sir_previous <- sampling(sir_mod, data = get_data_sir(nobs-1), iter = 2, chains = 1)
    fit_current  <- sampling(sir_mod, data = get_data_sir(nobs), iter = 2, chains = 1)
    initdist   <- get_sir_targetdist(fit_)
    targetdist <- get_sir_targetdist(fit_current)
    ## function to evaluate gradient of log density on unconstrained space
    results <- asmc_hmc(smctuning, targetdist, initdist, xparticles)
    print(results$infos_mcmc)
    xparticles <- results$particles$x
    xhistory[[length(xhistory)+1]] <- xparticles
    fit_  <- fit_current
  }
  xhistory
}
##
save(smctuning, smc_partial_results, nrep, file = "experiments/sirmodel/boardingschool-smc-partial.RData")
#

