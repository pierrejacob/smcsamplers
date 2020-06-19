rm(list = ls())
library(smcsamplers)
set.seed(1)
library(tidyverse)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
load("experiments/logistic/covtype.processed.RData")
#
dim(trainingset)

nrep <- 10
for (n in seq(from = 1e4, to = 1e5, by = 1e4)){
  p <- dim(trainingset)[2]-1
  Y <- trainingset[1:n,1] - 1
  X <- trainingset[1:n,2:(p+1)]
  ## standardization of inputs
  X <- apply(X, 2, function(v) (v - mean(v))/(sd(v)))
  X <- cbind(1, X)
  p <- p + 1

  ## Prior
  b <- matrix(0, nrow = p, ncol = 1)
  B <- diag(10, p, p)
  priordist <- get_mvnormal_diag(b[,1], diag(B))
  rinit <- function(n) b[,1] + sqrt(diag(B)) * matrix(rnorm(p*n), nrow = p)

  ## target distribution
  targetdist <- function(xs,...){
    pr <- priordist(xs)
    ll <- smcsamplers:::logistic_loglikelihood_gradient(xs, Y, X)
    return(list(logpdf = pr$logpdf + ll$logls, gradlogpdf = pr$gradlogpdf + ll$gradients))
  }
  ##
  ## find Laplace approximation
  objective <- function(par) -smcsamplers:::logistic_loglikelihood_gradient(matrix(par, ncol = 1), Y, X)$logls
  objective(rinit(1)[,1])
  optim_results <- optim(par = rinit(1)[,1], fn = objective, method = "BFGS")
  mle <- optim_results$par
  precision_mle <- numDeriv::hessian(objective, mle)
  variance_mle <- solve(precision_mle)
  ###

  ## SMC sampler from Laplace approximation to posterior
  laplacedist <- get_mvnormal(mle, variance_mle)
  ##
  smctuning <- list(nparticles = 2^10, ess_criterion = 0.5, nmoves = 1)
  smctuning$stepsize <- 0.3 * p^{-1/4}
  smctuning$nleapfrog <- ceiling(1/smctuning$stepsize)
  smc_results_laplace <- foreach(rep = 1:nrep) %dorng% {
    initparticles <- laplacedist$generate(smctuning$nparticles)
    asmc_hmc(smctuning, targetdist, laplacedist$eval, initparticles)
  }
  save(n, Y, X, b, B, smctuning, mle, precision_mle, variance_mle, smc_results_laplace,
       file = paste0("experiments/logistic/covtype.laplace.n", n, ".RData"))
}

# n <- 2e4
# load(file = paste0("experiments/logistic/covtype.laplace.n", n, ".RData"))
# adaptiveresults <- smc_results_laplace[[2]]
#
# nsteps <- length(adaptiveresults$xhistory)
# nsteps
#
#
# info_mcmc_df <- lapply(2:length(adaptiveresults$lambdas), function(ilambda) {
#   info_mcmc <- adaptiveresults$infos_mcmc[[ilambda-1]]
#   df_ <- data.frame(ilambda = ilambda, imove = 1:smctuning$nmoves,
#                     ar = sapply(info_mcmc, function(x) x$ar),
#                     sqjd = sapply(info_mcmc, function(x) x$sqjd))
# }) %>% bind_rows()
# ##
# gar <- ggplot(info_mcmc_df, aes(x = ilambda, y = ar)) + geom_point()
# gar <- gar + xlab("step") + ylab("acceptance rate") + geom_rangeframe()
# gsqjd <- ggplot(info_mcmc_df, aes(x = ilambda, y = sqjd)) + geom_point()
# gsqjd <- gsqjd + xlab("step") + ylab("relative jumping distance") + geom_rangeframe()
# # gridExtra::grid.arrange(gar, gsqjd, nrow = 2)
#
# sapply(adaptiveresults$roots_history, function(x) length(unique(x)))
# pathdf <- lapply(1:nsteps, function(time){
#   xs <- adaptiveresults$xhistory[[time]]
#   xmeans <- rowMeans(xs)
#   xvars <- apply(xs, 1, var)
#   data.frame(component = 1:p,
#              mean = xmeans,
#              var = xvars,
#              time = time)
# }) %>% bind_rows()
#
# ggplot(pathdf, aes(x = mean, y = var, group = interaction(component))) + geom_path() + scale_y_log10() +
#   geom_point(data=data.frame(mean = mle, var = diag(variance_mle)), aes(group=NULL)) +
#   geom_point(data = pathdf %>% filter(time == max(time)), aes(group = NULL), col = 'red')
#
# pathdf %>% filter(time == max(time)) %>% mutate(mle = mle)

