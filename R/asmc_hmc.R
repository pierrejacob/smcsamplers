## adaptive SMC sampler
## resample move when ESS is below a threshold
## HMC moves using diagonal mass matrix
## fitted on existing particles
#'@export
asmc_hmc <- function(smctuning, targetdist, initdist, initparticles){
  starttime <- Sys.time()
  particles <- list()
  particles$n <- smctuning$nparticles
  particles$x <- initparticles
  particles$d <- dim(particles$x)[1]
  particles$init <- initdist(particles$x)
  particles$target <- targetdist(particles$x)
  zerod <- rep(0, particles$d)
  Id <- diag(1, particles$d, particles$d)
  ### initialize the inverse temperature
  lambda_current <- 0
  lambdas <- c(lambda_current)
  ### initialize the log normalizing constant estimator
  log_ratio_normconst <- c()
  ess_realized <- c()
  ## history of particles
  xhistory <- list()
  istep <- 1
  xhistory[[istep]] <- particles$x
  ## root indices
  roots <- 1:particles$n
  roots_history <- list()
  roots_history[[1]] <- roots
  ahistory <- list()
  ahistory[[istep]] <- 1:particles$n
  # while inv temperature is not one...
  infos_mcmc <- list()
  while(lambda_current < 1){
    ## find the next inverse temperature
    ## function that maps lambda to a value in 1:N
    ess_deltalambda <- function(lambda){
      logw <- (lambda - lambda_current) * (particles$target$logpdf - particles$init$logpdf)
      logw[is.na(logw)] <- -Inf
      return(1/sum(smcsamplers::normalize_weight(logw)$nw^2))
    }
    ## if we can already set lambda = 1, do it
    if (ess_deltalambda(1) > smctuning$ess_criterion * smctuning$nparticles){
      lambda_next <- 1
    } else {
      ## otherwise use binary search
      search_results <- smcsamplers::search_lambda(lambda_current, ess_deltalambda, smctuning$ess_criterion * smctuning$nparticles)
      lambda_next <- search_results$x
    }
    ## compute log weights using next lambda
    logweights <- (lambda_next - lambda_current) * (particles$target$logpdf - particles$init$logpdf)
    logweights[is.na(logweights)] <- -Inf
    ## normalize weight in a numerically stable way
    nweights <- smcsamplers::normalize_weight(logweights)
    ## compute empirical variance of each component of the target
    estimated_means <- sapply(1:particles$d, function(component)
      sum(particles$x[component,] * nweights$nw))
    estimated_variances <- sapply(1:particles$d, function(component)
      sum((particles$x[component,]-estimated_means[component])^2 * nweights$nw))
    ## set mass "matrix" as inverse of variance (here only consider diagonal elements)
    mass_matrix <- 1/estimated_variances
    ## Cholesky decomposition
    mass_chol <- sqrt(mass_matrix)
    ## compute effective sample size
    ess_realized <- c(ess_realized, 1/sum(nweights$nw^2))
    ## compute normalizing constant ratio estimator
    log_ratio_normconst <- c(log_ratio_normconst, nweights$avew)
    ## multinomial resampling
    ancestors <- sample(x = 1:smctuning$nparticles, size = smctuning$nparticles,
                        prob = nweights$nw, replace = TRUE)
    ## keep track of roots of the genealogical tree
    roots <- roots[ancestors]
    ## perform resampling
    particles$x <- particles$x[,ancestors,drop=F]
    particles$init$logpdf <- particles$init$logpdf[ancestors]
    particles$init$gradlogpdf <- particles$init$gradlogpdf[,ancestors,drop=F]
    particles$target$logpdf <- particles$target$logpdf[ancestors]
    particles$target$gradlogpdf <- particles$target$gradlogpdf[,ancestors,drop=F]
    ## HMC move targeting distribution corresponding to next lambda
    info <- list()
    for (imove in 1:smctuning$nmoves){
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
    lambdas <- c(lambdas, lambda_current)
    istep <- istep + 1
    ## store history of particles and ancestry
    xhistory[[istep]] <- particles$x
    ahistory[[istep]] <- ancestors
    roots_history[[istep]] <- roots
  }
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(particles = particles, lambdas = lambdas, log_ratio_normconst = log_ratio_normconst,
              ess_realized = ess_realized, elapsedtime = elapsedtime, infos_mcmc = infos_mcmc,
              xhistory = xhistory, ahistory = ahistory, roots_history = roots_history))
}
