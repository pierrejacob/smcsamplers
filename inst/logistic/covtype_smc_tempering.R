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
#
n <- 1e3
ntest <- 1e3
p <- dim(trainingset)[2]-1
Ytest <- testset[1:ntest,1] - 1
Xtest <- testset[1:ntest,2:(p+1)]
Y <- trainingset[1:n,1] - 1
X <- trainingset[1:n,2:(p+1)]

summary(Y)
summary(X)

## standardization of inputs
X <- apply(X, 2, function(v) (v - mean(v))/(sd(v)))
X <- cbind(1, X)
Xtest <- apply(Xtest, 2, function(v) (v - mean(v))/(sd(v)))
Xtest <- cbind(1, Xtest)
p <- p + 1

dim(testset)

## Prior
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)

smctuning <- list(nparticles = 2^10, ess_criterion = 0.5, nmoves = 1)

## initial distribution
initdist <- get_mvnormal_diag(b[,1], diag(B))
initparticles <- b[,1] + sqrt(diag(B)) * matrix(rnorm(p*smctuning$nparticles), nrow = p)
## target distribution
targetdist <- function(xs,...){
  pr <- initdist(xs)
  ll <- smcsamplers:::logistic_loglikelihood_gradient(xs, Y, X)
  return(list(logpdf = pr$logpdf + ll$logls, gradlogpdf = pr$gradlogpdf + ll$gradients))
}
##

# targetdist(initparticles)

smctuning$stepsize <- 0.3 * p^{-1/4}
smctuning$nleapfrog <- ceiling(1/smctuning$stepsize)

adaptiveresults <- asmc_hmc(smctuning, targetdist, initdist, initparticles)
nsteps <- length(adaptiveresults$xhistory)

# qplot(x = 1:nsteps, y = adaptiveresults$lambdas, geom = "line") + ylab(expression(lambda)) + xlab("time")

## get information about how the MCMC steps performed
info_mcmc_df <- lapply(2:length(adaptiveresults$lambdas), function(ilambda) {
  info_mcmc <- adaptiveresults$infos_mcmc[[ilambda-1]]
  df_ <- data.frame(ilambda = ilambda, imove = 1:smctuning$nmoves,
                    ar = sapply(info_mcmc, function(x) x$ar),
                    sqjd = sapply(info_mcmc, function(x) x$sqjd))
}) %>% bind_rows()
##
gar <- ggplot(info_mcmc_df, aes(x = ilambda, y = ar)) + geom_point()
gar <- gar + xlab("step") + ylab("acceptance rate") + geom_rangeframe()
gsqjd <- ggplot(info_mcmc_df, aes(x = ilambda, y = sqjd)) + geom_point()
gsqjd <- gsqjd + xlab("step") + ylab("relative jumping distance") + geom_rangeframe()
# gridExtra::grid.arrange(gar, gsqjd, nrow = 2)
##

## finer grid of lambdas
nlambdas <- 50
require(cobs)
xgrid <- (0:(nsteps-1))/(nsteps-1)
fit = cobs(x = xgrid, y = adaptiveresults$lambdas, constraint= "increase",
           lambda=0, degree=1, # for L1 roughness
           knots=seq(0, 1, length.out=floor(nsteps/2)), # desired nr of knots
           tau=0.5) # to predict median
xfinergrid <- seq(from = 0, to = 1, length.out = nlambdas)
yrange <- predict(fit,interval="none",z=c(0,1))[,2]
yfinergrid <- (predict(fit,interval="none",z=xfinergrid)[,2] - yrange[1])/(yrange[2]-yrange[1])

# plot(xgrid, adaptiveresults$lambdas, xlim = c(0,1), ylim = c(0,1))
# points(xfinergrid, yfinergrid, col = 'red', lty = 2)
##
smctuning$lambdas <- yfinergrid
smctuning$lambdas[1] <- 0
smctuning$lambdas <- pmax(pmin(smctuning$lambdas, 1), 0)
##
# initparticles <- b[,1] + sqrt(diag(B)) * matrix(rnorm(p*smctuning$nparticles), nrow = p)
# results <- smc_hmc(smctuning, targetdist, initdist, initparticles)
# print(tail(results$nroots))
# qplot(x = 1:length(results$lambdas), results$nroots, geom = 'line') + xlab("time") + geom_rangeframe()


nrep <- 10
smc_hmc_results <- foreach(irep = 1:nrep) %dorng% {
  initparticles <- b[,1] + sqrt(diag(B)) * matrix(rnorm(p*smctuning$nparticles), nrow = p)
  smc_hmc(smctuning, targetdist, initdist, initparticles)
}
##
save(smctuning, nrep, smc_hmc_results, file = "experiments/logistic/covtype.smctempering.RData")




