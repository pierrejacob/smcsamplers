#'@export
smc_hmc <- function(smctuning, targetdist, initdist, initparticles){
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
  xmeans_history <- list()
  xvars_history <- list()
  xmeans_history[[1]] <- rowMeans(particles$x)
  xvars_history[[1]] <- apply(particles$x, 1, var)
  # while inv temperature is not one...
  while(lambda_current < 1){
    lambda_next <- lambdas[istep+1]
    ## now weight / resample / MALA move
    ## approximate variance of current target
    incrweight <- (lambda_next - lambda_current) * (particles$target$logpdf - particles$init$logpdf)
    incrweight[is.na(incrweight)] <- -Inf
    mlw <- max(incrweight)
    log_ratio_normconst <- c(log_ratio_normconst, mlw + log(sum(nweights * exp(incrweight - mlw))))
    logweights <- logweights + incrweight
    nweights <- PET::normalize_weight(logweights)$nw
    ess_realized <- c(ess_realized, 1/sum(nweights^2))
    ## compute empirical variance of each component of the target
    estimated_means <- sapply(1:particles$d, function(component)
      sum(particles$x[component,] * nweights))
    estimated_variances <- sapply(1:particles$d, function(component)
      sum((particles$x[component,]-estimated_means[component])^2 * nweights))
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
    xmeans_history[[istep]] <- estimated_means
    xvars_history[[istep]] <- estimated_variances
    nroots <- c(nroots, length(unique(roots)))
  }
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(particles = particles, lambdas = lambdas,
              log_ratio_normconst = log_ratio_normconst, ess_realized = ess_realized,
              elapsedtime = elapsedtime, xmeans_history = xmeans_history, xvars_history = xvars_history,
              nroots = nroots, infos_mcmc = infos_mcmc)) #xhistory = xhistory, nwhistory = nwhistory,
  # ahistory = ahistory, roots_history = roots_history
}
