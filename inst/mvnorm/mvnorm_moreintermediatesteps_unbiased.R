#
rm(list = ls())
library(smcsamplers)
graphsettings <- set_custom_theme()
registerDoParallel(cores = 8)
set.seed(3)
#
smc_hmc_nonadaptive <- function(smctuning, targetdist, initdist, initparticles){
  starttime <- Sys.time()
  particles <- list()
  particles$x <- initparticles
  particles$init <- initdist(particles$x)
  particles$target <- targetdist(particles$x)
  particles$n <- smctuning$nparticles
  particles$d <- dim(particles$x)[1]
  zerod <- rep(0, particles$d)
  Id <- diag(1, particles$d, particles$d)
  ### initialize the inverse temperature
  lambdas <- smctuning$lambdas
  lambda_current <- lambdas[1]
  ## lambda_current should be = 0
  ## initialize the log normalizing constant estimator
  log_ratio_normconst <- c()
  ##
  ess_realized <- c()
  istep <- 1
  roots <- 1:particles$n
  logweights <- rep(0, particles$n)
  nweights <- rep(1/particles$n, particles$n)
  infos_mcmc <- list()
  nroots <- c(particles$n)
  moment1_history <- list()
  moment2_history <- list()
  moment1_history[[1]] <- rowMeans(particles$x)
  moment2_history[[1]] <- (particles$x %*% t(particles$x)) / smctuning$nparticles
  # while inv temperature is not one...
  while(lambda_current < 1){
    lambda_next <- lambdas[istep+1]
    ## now weight / resample / HMC move
    ## approximate variance of current target
    incrweight <- (lambda_next - lambda_current) * (particles$target$logpdf - particles$init$logpdf)
    incrweight[is.na(incrweight)] <- -Inf
    mlw <- max(incrweight)
    log_ratio_normconst <- c(log_ratio_normconst, mlw + log(sum(nweights * exp(incrweight - mlw))))
    logweights <- logweights + incrweight
    nweights <- smcsamplers::normalize_weight(logweights)$nw
    ess_realized <- c(ess_realized, 1/sum(nweights^2))
    estimated_variances <- smctuning$variances[[istep+1]]
    ## set mass "matrix" as inverse of variance (here only consider diagonal elements)
    mass_matrix <- 1/estimated_variances
    ## Cholesky decomposition
    mass_chol <- sqrt(mass_matrix)

    if (1/sum(nweights^2) < (smctuning$ess_criterion * smctuning$nparticles)){
      ## resampling
      ancestors <- sample(x = 1:smctuning$nparticles, size = smctuning$nparticles, prob = nweights, replace = TRUE)
      logweights <- rep(0, particles$n)
      nweights <- rep(1/particles$n, particles$n)
      roots <- roots[ancestors]
      particles$x <- particles$x[,ancestors,drop=F]
      particles$init$logpdf <- particles$init$logpdf[ancestors]
      particles$init$gradlogpdf <- particles$init$gradlogpdf[,ancestors,drop=F]
      particles$target$logpdf <- particles$target$logpdf[ancestors]
      particles$target$gradlogpdf <- particles$target$gradlogpdf[,ancestors,drop=F]
    } else {
      ancestors <- 1:smctuning$nparticles
    }
    info <- list()
    for (imove in 1:smctuning$nmoves){
      ## HMC move
      lambda <- lambda_next
      ## draw momenta variables
      initial_momenta <- matrix(rnorm(particles$n * particles$d), nrow = particles$d) * mass_chol
      ## compute gradients wrt target at current positions
      grads_ <- (1-lambda) * particles$init$gradlogpdf + lambda * particles$target$gradlogpdf
      positions <- particles$x
      ## leap frog integrator
      momenta <- initial_momenta + smctuning$stepsize * grads_ / 2
      for (step in 1:smctuning$nleapfrog){
        positions <- positions + smctuning$stepsize * momenta / mass_matrix
        eval_init <- initdist(positions)
        eval_target <- targetdist(positions)
        if (step != smctuning$nleapfrog){
          momenta <- momenta + smctuning$stepsize * ((1-lambda) * eval_init$gradlogpdf + lambda * eval_target$gradlogpdf)
        }
      }
      momenta <- momenta + smctuning$stepsize * ((1-lambda) * eval_init$gradlogpdf + lambda * eval_target$gradlogpdf) / 2
      ## compute MH acceptance ratios
      proposed_pdfs <- (1-lambda) * eval_init$logpdf + lambda * eval_target$logpdf
      current_pdfs  <- (1-lambda) * particles$init$logpdf + lambda * particles$target$logpdf
      mhratios <- proposed_pdfs - current_pdfs
      mhratios <- mhratios + (-0.5*colSums((momenta/mass_chol)^2)) - (-0.5*colSums((initial_momenta/mass_chol)^2))
      if (any(is.na(mhratios))) mhratios[is.na(mhratios)] <- -Inf
      ## compute expected squared jumping distance
      ## normalized by estimated variance component-wise
      sqjd <- mean(pmin(1, exp(mhratios)) * colSums((particles$x - positions)^2/estimated_variances)/particles$d)
      ## accept - reject step
      accepts <- log(runif(particles$n)) < mhratios
      ## replace accepted particles
      particles$x[,accepts] <- positions[,accepts,drop=F]
      particles$init$logpdf[accepts] <- eval_init$logpdf[accepts]
      particles$target$logpdf[accepts] <- eval_target$logpdf[accepts]
      particles$target$gradlogpdf[,accepts] <- eval_target$gradlogpdf[,accepts,drop=F]
      particles$init$gradlogpdf[,accepts] <- eval_init$gradlogpdf[,accepts,drop=F]
      ## report performance of MCMC moves
      info[[imove]] <- list(ar = mean(accepts), sqjd = sqjd)
    }
    infos_mcmc[[istep]] <- info
    lambda_current <- lambda_next
    istep <- istep + 1
    ## store mean and variance
    ## compute empirical variance of each component of the target
    estimated_moment1 <- sapply(1:particles$d, function(component)
      sum(particles$x[component,] * nweights))
    # estimated_moment2 <- sapply(1:particles$d, function(component)
    # sum(particles$x[component,]^2 * nweights))
    estimated_moment2 <- (particles$x %*% (nweights * t(particles$x)))
    moment1_history[[istep]] <- estimated_moment1
    moment2_history[[istep]] <- estimated_moment2
    nroots <- c(nroots, length(unique(roots)))
  }
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(particles = particles, lambdas = lambdas,
              log_ratio_normconst = log_ratio_normconst, ess_realized = ess_realized,
              elapsedtime = elapsedtime, moment1_history = moment1_history, moment2_history = moment2_history,
              nroots = nroots, infos_mcmc = infos_mcmc))
}

## set up problem
nrep <- 100
dimensions <- c(32, 64, 128, 256, 512)
# dimensions <- c(32, 64, 128)
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
    sqerror_means <- mean((results$moment1_history[[nsteps]] - means[,nsteps])^2)
    zhat <- sum(results$log_ratio_normconst)
    nroots <- results$nroots[nsteps]
    return(c(nsteps, sqerror_means, zhat, nroots))
  }

  indep_results <- foreach(irep = 1:nrep) %dorng% {
    initparticles <- initmean + sqrt(initvar) * matrix(rnorm(dimension*smctuning$nparticles), nrow = dimension)
    adaptive_smc_run <- smc_hmc(smctuning, targetdist, initdist, initparticles)
    smctuning$variances <- adaptive_smc_run$xvars_history
    nonadaptive_smc_run <- smc_hmc_nonadaptive(smctuning, targetdist, initdist, initparticles)
    process_smc_output(nonadaptive_smc_run)
  }
  indep_results_df <- data.frame(do.call(what = rbind, args = indep_results))
  names(indep_results_df) <- c("nsteps", "sqerror", "zhat", "nroots")
  indep_results_df$rep <- 1:nrep
  indep_results_df$nparticles <- smctuning$nparticles
  indep_results_df$dimension <- dimension
  indep_results_df
  intermedf <- rbind(intermedf, indep_results_df)
}

intermedf.unbiased <- intermedf
save(intermedf.unbiased, dimensions, smctuning, file = "experiments/mvnorm/scalingintermediate.unbiased.RData")

