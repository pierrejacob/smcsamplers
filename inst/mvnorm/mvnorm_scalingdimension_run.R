## script to study impact of dimension
## initial distribution and target distribution are both multivariate Normals
## geometric bridge between these distributions
## SMC samplers is run with adaptive selection of inverse temperature lambda and resampling
## HMC moves are performed with adaptation of the mass matrix

rm(list = ls())
library(smcsamplers)
graphsettings <- set_custom_theme()
library(doParallel)
library(doRNG)
library(tidyverse)
registerDoParallel(cores = 6)
set.seed(3)

smctuning <- list()
smctuning$ess_criterion <- 0.5
smctuning$nmoves <- 1
smctuning$nparticles <- 2^8

nrep <- 50
dimensions <- c(32, 64, 128, 256, 512)
fixedndf <- data.frame()
for (idim in seq_along(dimensions)){
  dimension <- dimensions[idim]
  cat("dimension", dimension, "\n")
  initmean <- rep(1, dimension)
  initvar <- rep(0.5, dimension)
  targetmean <- rep(0, dimension)
  targetvar <- rep(1, dimension)
  initdist <- get_mvnormal_diag(initmean, initvar)
  targetdist <- get_mvnormal_diag(targetmean, targetvar)
  ## get exact mean and variance of intermediate distributions
  meanvar_intermediate <- function(lambda){
    precision <- (1-lambda) * (1/initvar) + lambda * (1/targetvar)
    var <- 1/precision
    mean <- var * ((1-lambda) * (1/initvar) * initmean + lambda * (1/targetvar) * targetmean)
    return(list(mean = mean, var = var))
  }
  ##
  ## process output of smc algorithm
  process_smc_output <- function(results){
    nsteps <- length(results$xhistory)
    means <- sapply(results$lambdas, function(lambda) meanvar_intermediate(lambda)$mean)
    sqerror_means <- mean((rowMeans(results$particles$x) - means[,nsteps])^2)
    zhat <- sum(results$log_ratio_normconst)
    nroots <- length(unique(results$roots_history[[nsteps]]))
    return(c(nsteps, sqerror_means, zhat, nroots))
  }
  ##
  smctuning$stepsize <- 1 * dimension^{-1/4}
  smctuning$nleapfrog <- ceiling(1/smctuning$stepsize)
  ##
  indep_results <- foreach(irep = 1:nrep) %dorng% {
    initparticles <- initmean + sqrt(initvar) * matrix(rnorm(dimension*smctuning$nparticles), nrow = dimension)
    process_smc_output(asmc_hmc(smctuning, targetdist, initdist, initparticles))
  }
  indep_results_df <- data.frame(do.call(what = rbind, args = indep_results))
  names(indep_results_df) <- c("nsteps", "sqerror", "zhat", "nroots")
  indep_results_df$rep <- 1:nrep
  indep_results_df$nparticles <- smctuning$nparticles
  indep_results_df$dimension <- dimension
  fixedndf <- rbind(fixedndf, indep_results_df)
}

save(fixedndf, dimensions, smctuning, file = "experiments/mvnorm/scalingfixedn.RData")

# head(fixedndf)
# gnsteps <- ggplot(fixedndf, aes(x = dimension, y = nsteps, group = rep)) + geom_point()
# gsqerror <- ggplot(fixedndf  %>% group_by(dimension) %>% summarize(meansqerror = mean(sqerror)),
#                    aes(x = dimension, y = meansqerror)) + geom_point()
# gzhat <- ggplot(fixedndf %>% group_by(dimension) %>% summarise(varzhat = var(zhat)),
#                 aes(x = dimension, y = varzhat)) + geom_point()
# gnroots <- ggplot(fixedndf, aes(x = dimension, y = nroots, group = rep)) + geom_point()
#
# print(gridExtra::grid.arrange(gnsteps, gnroots, gsqerror, gzhat))

linearndf <- data.frame()
for (idim in seq_along(dimensions)){
  dimension <- dimensions[idim]
  cat("dimension", dimension, "\n")
  initmean <- rep(1, dimension)
  initvar <- rep(0.5, dimension)
  targetmean <- rep(0, dimension)
  targetvar <- rep(1, dimension)
  initdist <- get_mvnormal_diag(initmean, initvar)
  targetdist <- get_mvnormal_diag(targetmean, targetvar)
  ## get exact mean and variance of intermediate distributions
  meanvar_intermediate <- function(lambda){
    precision <- (1-lambda) * (1/initvar) + lambda * (1/targetvar)
    var <- 1/precision
    mean <- var * ((1-lambda) * (1/initvar) * initmean + lambda * (1/targetvar) * targetmean)
    return(list(mean = mean, var = var))
  }
  ##
  ## process output of smc algorithm
  process_smc_output <- function(results){
    nsteps <- length(results$xhistory)
    means <- sapply(results$lambdas, function(lambda) meanvar_intermediate(lambda)$mean)
    sqerror_means <- mean((rowMeans(results$particles$x) - means[,nsteps])^2)
    zhat <- sum(results$log_ratio_normconst)
    nroots <- length(unique(results$roots_history[[nsteps]]))
    return(c(nsteps, sqerror_means, zhat, nroots))
  }
  smctuning$stepsize <- 1 * dimension^{-1/4}
  smctuning$nleapfrog <- ceiling(1/smctuning$stepsize)
  smctuning$nparticles <- 2^8 + 8 * dimension
  ##
  indep_results <- foreach(irep = 1:nrep) %dorng% {
    initparticles <- initmean + sqrt(initvar) * matrix(rnorm(dimension*smctuning$nparticles), nrow = dimension)
    process_smc_output(asmc_hmc(smctuning, targetdist, initdist, initparticles))
  }
  indep_results_df <- data.frame(do.call(what = rbind, args = indep_results))
  names(indep_results_df) <- c("nsteps", "sqerror", "zhat", "nroots")
  indep_results_df$rep <- 1:nrep
  indep_results_df$nparticles <- smctuning$nparticles
  indep_results_df$dimension <- dimension
  linearndf <- rbind(linearndf, indep_results_df)
}
save(linearndf, dimensions, smctuning, file = "experiments/mvnorm/scalinglinearn.RData")

# gnsteps <- ggplot(linearndf, aes(x = dimension, y = nsteps, group = rep)) + geom_point()
# gsqerror <- ggplot(linearndf  %>% group_by(dimension) %>% summarize(meansqerror = mean(sqerror)),
#                    aes(x = dimension, y = meansqerror)) + geom_point()
# gzhat <- ggplot(linearndf %>% group_by(dimension) %>% summarise(varzhat = var(zhat)),
#                 aes(x = dimension, y = varzhat)) + geom_point()
# gnroots <- ggplot(linearndf, aes(x = dimension, y = nroots, group = rep)) + geom_point()
#
# print(gridExtra::grid.arrange(gnsteps, gnroots, gsqerror, gzhat))
