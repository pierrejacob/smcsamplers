## smctuning list containing
## targetdist
## initdist
## initparticles
#'@export
asmc_mala <- function(smctuning, targetdist, initdist, initparticles){
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
  lambda_current <- 0
  lambdas <- c(lambda_current)
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
  # while inv temperature is not one...
  while(lambda_current < 1){
    ## find the next inverse temperature
    ## based on standard IS
    ess_deltalambda <- function(lambda){
      logw <- (lambda - lambda_current) * (particles$target$logpdf - particles$init$logpdf)
      logw[is.na(logw)] <- -Inf
      return(1/sum(smcsamplers::normalize_weight(logw)$nw^2))
    }
    if (ess_deltalambda(1) > smctuning$ess_criterion * smctuning$nparticles){
      lambda_next <- 1
    } else {
      search_results <- smcsamplers::search_lambda(lambda_current, ess_deltalambda, smctuning$ess_criterion * smctuning$nparticles)
      lambda_next <- search_results$x
    }
    ## now weight / resample / MALA move
    ## approximate variance of current target
    logweights <- (lambda_next - lambda_current) * (particles$target$logpdf - particles$init$logpdf)
    logweights[is.na(logweights)] <- -Inf
    nweights <- smcsamplers::normalize_weight(logweights)
    plainIS_weightedcov <- cov.wt(x = t(particles$x), wt = nweights$nw)$cov
    smctuning$cov <- plainIS_weightedcov
    smctuning$cholcov <- chol(smctuning$cov)
    smctuning$cholinvcov <- t(chol(solve(smctuning$cov)))
    ess_realized <- c(ess_realized, 1/sum(nweights$nw^2))
    log_ratio_normconst <- c(log_ratio_normconst, nweights$avew)
    ## resampling
    ancestors <- sample(x = 1:smctuning$nparticles, size = smctuning$nparticles, prob = nweights$nw, replace = TRUE)
    eves <- eves[ancestors]
    particles$x <- particles$x[,ancestors,drop=F]
    particles$init$logpdf <- particles$init$logpdf[ancestors]
    particles$init$gradlogpdf <- particles$init$gradlogpdf[,ancestors,drop=F]
    particles$target$logpdf <- particles$target$logpdf[ancestors]
    particles$target$gradlogpdf <- particles$target$gradlogpdf[,ancestors,drop=F]
    ## MALA move
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
      MHratios <- (1-lambda) * next_particles$init$logpdf + lambda * next_particles$target$logpdf
      MHratios <- MHratios - ((1-lambda) * particles$init$logpdf + lambda * particles$target$logpdf)
      MHratios <- MHratios + bw - fw
      MHratios[is.na(MHratios)] <- -Inf
      for (iparticle in 1:particles$n){
        if (log(us[iparticle]) < MHratios[iparticle]){
          particles$x[,iparticle] <- next_particles$x[,iparticle]
          particles$init$logpdf[iparticle] <- next_particles$init$logpdf[iparticle]
          particles$init$gradlogpdf[,iparticle] <- next_particles$init$gradlogpdf[,iparticle]
          particles$target$logpdf[iparticle] <- next_particles$target$logpdf[iparticle]
          particles$target$gradlogpdf[,iparticle] <- next_particles$target$gradlogpdf[,iparticle]
        }
      }
    }
    lambda_current <- lambda_next
    lambdas <- c(lambdas, lambda_current)
    istep <- istep + 1
    xhistory[[istep]] <- particles$x
    ahistory[[istep]] <- ancestors
    eves_history[[istep]] <- eves
  }
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(particles = particles, lambdas = lambdas,
              log_ratio_normconst = log_ratio_normconst, ess_realized = ess_realized,
              elapsedtime = elapsedtime, xhistory = xhistory, ahistory = ahistory, eves_history = eves_history))
}
