## adaptive SMCS using MALA moves
## resampling is done according to ESS criterion
## sequence of lambdas is fixed,
## and MCMC moves are fixed too
## smctuning list containing
## - lambdas, nparticles, ess_criterion, nmoves, sqrth, h
## targetdist: describes the target distribution
## initdist: describes the initial distribution
## initparticles:
#'@export
smc_mala <- function(smctuning, targetdist, initdist, initparticles){
  starttime <- Sys.time()
  particles <- list()
  particles$x <- initparticles
  particles$init <- initdist(particles$x)
  particles$target <- targetdist(particles$x)
  particles$n <- smctuning$nparticles
  particles$d <- dim(particles$x)[1]
  zerod <- rep(0, particles$d)
  Id <- diag(1, particles$d, particles$d)
  ##
  ### initialize the inverse temperature
  lambdas <- smctuning$lambdas
  lambda_current <- lambdas[1]
  ## lambda_current should be = 0
  ### initialize the log normalizing constant estimator
  log_ratio_normconst <- c()
  ##
  ess_realized <- c()
  ## history of particles
  xhistory <- list()
  istep <- 1
  xhistory[[istep]] <- particles$x
  ## eve indices
  eves <- 1:particles$n
  eves_history <- list()
  eves_history[[1]] <- eves
  ahistory <- list()
  ahistory[[istep]] <- 1:particles$n
  logweights <- rep(0, particles$n)
  nweights <- rep(1/particles$n, particles$n)
  nwhistory <- list()
  nwhistory[[1]] <- nweights
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
    nweights <- smcsamplers::normalize_weight(logweights)$nw
    ess_realized <- c(ess_realized, 1/sum(nweights^2))
    ## compute covariance of target
    smctuning$cov <- cov.wt(x = t(particles$x), wt = nweights)$cov
    smctuning$cholcov <- chol(smctuning$cov)
    smctuning$cholinvcov <- t(chol(solve(smctuning$cov)))
    if (1/sum(nweights^2) < (smctuning$ess_criterion * smctuning$nparticles)){
      ## resampling
      ancestors <- sample(x = 1:smctuning$nparticles, size = smctuning$nparticles, prob = nweights, replace = TRUE)
      logweights <- rep(0, particles$n)
      nweights <- rep(1/particles$n, particles$n)
      eves <- eves[ancestors]
      particles$x <- particles$x[,ancestors,drop=F]
      particles$init$logpdf <- particles$init$logpdf[ancestors]
      particles$init$gradlogpdf <- particles$init$gradlogpdf[,ancestors,drop=F]
      particles$target$logpdf <- particles$target$logpdf[ancestors]
      particles$target$gradlogpdf <- particles$target$gradlogpdf[,ancestors,drop=F]
    } else {
      ancestors <- 1:smctuning$nparticles
    }
    for (imove in 1:smctuning$nmoves){
      ## MALA move
      lambda <- lambda_next
      noise_ <- smctuning$sqrth * t(smcsamplers:::rmvnorm_cholesky_(particles$n, zerod, smctuning$cholcov))
      grads_ <- (1-lambda) * particles$init$gradlogpdf + lambda * particles$target$gradlogpdf
      ## X' = X + 0.5 * h Sigma %*% grad(X) + Normal(0, h Sigma)
      next_particles <- list(x = particles$x + 0.5 * smctuning$h * (smctuning$cov %*% grads_) + noise_,
                             n = particles$n, d = particles$d)
      ## target and initial densities at new points X'
      next_particles$target <- targetdist(next_particles$x)
      next_particles$init <- initdist(next_particles$x)
      ## forward weight: Normal(X'; X + 0.5 * h Sigma %*% grad(X), h Sigma)
      fw <- smcsamplers:::dmvnorm_cholesky_inverse(t(noise_), mean = zerod, cholesky_inverse = 1/smctuning$sqrth * smctuning$cholinvcov)
      ## backward weight: Normal(X; X' + 0.5 * h Sigma %*% grad(X'), h Sigma)
      back_grads_ <- (1-lambda) * next_particles$init$gradlogpdf + lambda * next_particles$target$gradlogpdf
      bw <- smcsamplers:::dmvnorm_cholesky_inverse(t(particles$x-(next_particles$x + 0.5 * smctuning$h * (smctuning$cov %*%  back_grads_))),
                                                   mean = zerod, cholesky_inverse = 1/smctuning$sqrth * smctuning$cholinvcov)
      ## accept-reject step
      us <- runif(particles$n)
      MRTHratios <- (1-lambda) * next_particles$init$logpdf + lambda * next_particles$target$logpdf
      MRTHratios <- MRTHratios - ((1-lambda) * particles$init$logpdf + lambda * particles$target$logpdf)
      MRTHratios <- MRTHratios + bw - fw
      MRTHratios[is.na(MRTHratios)] <- -Inf
      for (iparticle in 1:particles$n){
        if (log(us[iparticle]) < MRTHratios[iparticle]){
          particles$x[,iparticle] <- next_particles$x[,iparticle]
          particles$init$logpdf[iparticle] <- next_particles$init$logpdf[iparticle]
          particles$init$gradlogpdf[,iparticle] <- next_particles$init$gradlogpdf[,iparticle]
          particles$target$logpdf[iparticle] <- next_particles$target$logpdf[iparticle]
          particles$target$gradlogpdf[,iparticle] <- next_particles$target$gradlogpdf[,iparticle]
        }
      }
    }
    lambda_current <- lambda_next
    istep <- istep + 1
    xhistory[[istep]] <- particles$x
    ahistory[[istep]] <- ancestors
    eves_history[[istep]] <- eves
    nwhistory[[istep]] <- nweights
  }
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(particles = particles, lambdas = lambdas,
              log_ratio_normconst = log_ratio_normconst, ess_realized = ess_realized,
              elapsedtime = elapsedtime, xhistory = xhistory, nwhistory = nwhistory, ahistory = ahistory, eves_history = eves_history))
}
