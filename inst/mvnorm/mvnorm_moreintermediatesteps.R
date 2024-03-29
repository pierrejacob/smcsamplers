## script to see what happens if one adds more and more intermediate distributions
## initial distribution and target distribution are both multivariate Normals
## geometric bridge between these distributions
## SMC samplers is run with adaptive resampling
## HMC moves are performed with adaptation of the mass matrix

rm(list = ls())
library(smcsamplers)
graphsettings <- set_custom_theme()
registerDoParallel(cores = 8)
set.seed(3)

## set up problem
nrep <- 100
dimensions <- c(32, 64, 128, 256, 512)
intermedf <- data.frame()
for (idim in seq_along(dimensions)){
  dimension <- dimensions[idim]
  cat("dimension", dimension, "\n")
  nintermediate <- dimension * 1
  initmean <- rep(1, dimension)
  initvar <- rep(0.5, dimension)
  targetmean <- rep(0, dimension)
  targetvar <- rep(1, dimension)
  initdist <- get_mvnormal_diag(initmean, initvar)
  targetdist <- get_mvnormal_diag(targetmean, targetvar)
  ## get exact mean and variance of intermediate distribution...
  meanvar_intermediate <- function(lambda){
    precision <- (1-lambda) * (1/initvar) + lambda * (1/targetvar)
    var <- 1/precision
    mean <- var * ((1-lambda) * (1/initvar) * initmean + lambda * (1/targetvar) * targetmean)
    return(list(mean = mean, var = var))
  }
  ## set tuning parameters for SMC sampler
  smctuning <- list()
  ## number of particles
  smctuning$nparticles <- 2^8
  ## ESS criterion (number in (0,1))
  smctuning$ess_criterion <- 0.5
  ## stepsize for the HMC moves
  smctuning$stepsize <- 1 * dimension^{-1/4}
  ## number of leapfrog steps
  smctuning$nleapfrog <- ceiling(1/smctuning$stepsize)
  ## number of HMC moves to perform
  smctuning$nmoves <- 1
  ## one run:
  initparticles <- initmean + sqrt(initvar) * matrix(rnorm(dimension*smctuning$nparticles), nrow = dimension)
  results <- asmc_hmc(smctuning, targetdist, initdist, initparticles)
  ## what's in the results?
  nsteps <- length(results$xhistory)

  ## finer grid of lambdas
  nlambdas <- nintermediate
  require(cobs)
  xgrid <- (0:(nsteps-1))/(nsteps-1)
  fit = cobs(x = xgrid, y = results$lambdas, constraint= "increase",
             lambda=0, degree=1, # for L1 roughness
             knots=seq(0, 1, length.out=floor(nsteps/2)), # desired nr of knots
             tau=0.5) # to predict median
  xfinergrid <- seq(from = 0, to = 1, length.out = nlambdas)
  yrange <- predict(fit,interval="none",z=c(0,1))[,2]
  yfinergrid <- (predict(fit,interval="none",z=xfinergrid)[,2] - yrange[1])/(yrange[2]-yrange[1])

  # plot(xgrid, results$lambdas, xlim = c(0,1), ylim = c(0,1))
  # points(xfinergrid, yfinergrid, col = 'red', lty = 2)

  ## set up grid of lambdas
  smctuning$lambdas <- yfinergrid
  smctuning$lambdas[1] <- 0
  smctuning$lambdas <- pmax(pmin(smctuning$lambdas, 1), 0)

  ## process output of smc algorithm
  process_smc_output <- function(results){
    nsteps <- length(results$lambdas)
    means <- sapply(results$lambdas, function(lambda) meanvar_intermediate(lambda)$mean)
    sqerror_means <- mean((results$xmeans_history[[nsteps]] - means[,nsteps])^2)
    zhat <- sum(results$log_ratio_normconst)
    nroots <- results$nroots[nsteps]
    return(c(nsteps, sqerror_means, zhat, nroots))
  }

  indep_results <- foreach(irep = 1:nrep) %dorng% {
    initparticles <- initmean + sqrt(initvar) * matrix(rnorm(dimension*smctuning$nparticles), nrow = dimension)
    process_smc_output(smc_hmc(smctuning, targetdist, initdist, initparticles))
  }
  indep_results_df <- data.frame(do.call(what = rbind, args = indep_results))
  names(indep_results_df) <- c("nsteps", "sqerror", "zhat", "nroots")
  indep_results_df$rep <- 1:nrep
  indep_results_df$nparticles <- smctuning$nparticles
  indep_results_df$dimension <- dimension
  indep_results_df
  intermedf <- rbind(intermedf, indep_results_df)
}

save(intermedf, dimensions, smctuning, file = "experiments/mvnorm/scalingintermediate.RData")

# load("experiments/mvnorm/scalingintermediate.RData")
# tail(intermedf)
# load("experiments/mvnorm/scalingintermediate.unbiased.RData")
# tail(intermedf.unbiased)
#
# df_ <- rbind(intermedf %>% mutate(case = "fixed N & d steps"),
#              intermedf.unbiased %>% mutate(case = "fixed N & d steps - unbiased"))
# head(df_)
#
# plotwidth <- 7
# plotheight <- 5
# linesize <- 1
#
# colours <- c("antiquewhite3", "orange")
#
# gnroots <- ggplot(df_ %>% group_by(dimension, case) %>% summarise(meannroots = mean(nroots)),
#                   aes(x = dimension, y = meannroots, colour = case, linetype = case)) + geom_line(size = linesize)
# gnroots <- gnroots + scale_color_manual(name = '', values = colours[c(1,2)]) + ylab("mean # roots")
# gnroots <- gnroots + scale_linetype(name = "")
# gnroots <- gnroots + geom_rangeframe() + scale_x_continuous(breaks = dimensions)
# gnroots
#
#
# gsqerror <- ggplot(df_  %>% group_by(dimension, case) %>% summarize(meansqerror = mean(sqerror)),
#                    aes(x = dimension, y = meansqerror, colour = case, linetype = case)) + geom_line(size = linesize)
# gsqerror <- gsqerror + scale_color_manual(name = '', values = colours[c(1,2)]) + ylab("MSE E[X]") + scale_linetype(name = "")
# gsqerror <- gsqerror + geom_rangeframe()  + scale_x_continuous(breaks = dimensions)
# gsqerror
#
# gzhat <- ggplot(df_ %>% group_by(dimension, case) %>% summarise(varzhat = var(zhat)),
#                 aes(x = dimension, y = varzhat, colour = case, linetype = case)) + geom_line(size = linesize)
# gzhat <- gzhat + scale_color_manual(name = '', values = colours[c(1,2)]) + ylab("variance log constant")  + scale_linetype(name = "")
# gzhat <- gzhat + geom_rangeframe() + scale_x_continuous(breaks = dimensions)
# gzhat
#
#
# gmsezhat <- ggplot(df_ %>% group_by(dimension, case) %>% summarise(mse_logzhat = mean((exp(zhat)-1)^2)),
#                    aes(x = dimension, y = mse_logzhat, colour = case, linetype = case)) + geom_line(size = linesize)
# gmsezhat <- gmsezhat + scale_color_manual(name = '', values = colours[c(1,2)]) + ylab("MSE constant")  + scale_linetype(name = "")
# gmsezhat <- gmsezhat + geom_rangeframe() + scale_x_continuous(breaks = dimensions)
# gmsezhat # + scale_y_log10()
#
#
# gmsezhat <- ggplot(df_ %>% group_by(dimension, case) %>% summarise(mse_logzhat = mean((zhat-0)^2)),
#                    aes(x = dimension, y = mse_logzhat, colour = case, linetype = case)) + geom_line(size = linesize)
# gmsezhat <- gmsezhat + scale_color_manual(name = '', values = colours[c(1,2)]) + ylab("MSE log constant")  + scale_linetype(name = "")
# gmsezhat <- gmsezhat + geom_rangeframe() + scale_x_continuous(breaks = dimensions)
# gmsezhat # + scale_y_log10()

