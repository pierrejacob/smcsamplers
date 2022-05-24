rm(list = ls())
library(smcsamplers)
set.seed(1)
library(tidyverse)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
library(ggridges)
library(ggthemes)
theme_set(theme_tufte(ticks = TRUE))
theme_update(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
             axis.title.x = element_text(size = 25, margin = margin(20, 0, 0, 0), hjust = 1),
             axis.title.y = element_text(size = 25, angle = 90, margin = margin(0, 20, 0, 0), vjust = 1),
             legend.text = element_text(size = 20),
             legend.title = element_text(size = 20), title = element_text(size = 30),
             strip.text = element_text(size = 25), strip.background = element_rect(fill = "white"),
             legend.position = "bottom")

load("experiments/logistic/covtype.processed.RData")
trainingset <- rbind(trainingset, testset)
## SMC with HMC moves
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
    ## now weight / resample / MALA move
    ## approximate variance of current target
    incrweight <- (lambda_next - lambda_current) * (particles$target$logpdf - particles$init$logpdf)
    incrweight[is.na(incrweight)] <- -Inf
    mlw <- max(incrweight)
    log_ratio_normconst <- c(log_ratio_normconst, mlw + log(sum(nweights * exp(incrweight - mlw))))
    logweights <- logweights + incrweight
    nweights <- PET::normalize_weight(logweights)$nw
    ess_realized <- c(ess_realized, 1/sum(nweights^2))
    estimated_variances <- diag(smctuning$variances[[istep+1]])
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
              nroots = nroots, infos_mcmc = infos_mcmc)) #xhistory = xhistory, nwhistory = nwhistory,
  # ahistory = ahistory, roots_history = roots_history
}

unbiased_moments <- function(initdist, targetdist, smctuning){
  ## run SMC
  # initparticles1 <- initmean + sqrt(initvar) * matrix(rnorm(p*smctuning$nparticles), nrow = p)
  initparticles1 <- initdist$generate(smctuning$nparticles)
  smc_results1 <- smc_hmc_nonadaptive(smctuning, targetdist, initdist$eval, initparticles1)
  nsteps1 <- length(smc_results1$lambdas)
  moment11 <- smc_results1$moment1_history[[nsteps1]]
  moment21  <- smc_results1$moment2_history[[nsteps1]]
  logzhat1 <- sum(smc_results1$log_ratio_normconst)
  unbiased_moment1 <- moment11
  unbiased_moment2 <- moment21
  ## run second SMC
  initparticles2 <- initdist$generate(smctuning$nparticles)
  smc_results2 <- smc_hmc_nonadaptive(smctuning, targetdist, initdist$eval, initparticles2)
  nsteps2 <- length(smc_results2$lambdas)
  moment12 <- smc_results2$moment1_history[[nsteps2]]
  moment22  <- smc_results2$moment2_history[[nsteps2]]
  logzhat2 <- sum(smc_results2$log_ratio_normconst)
  ## propose second to first chain
  meetingtime <- Inf
  iteration <- 1
  meeting <- FALSE
  if (log(runif(1)) < (logzhat2 - logzhat1)){
    ## accept
    meeting <- TRUE
    meetingtime <- iteration
    logzhat1 <- logzhat2
    ## we can return an unbiased estimator
  } else {
    ## add bias correction term
    unbiased_moment1 <- unbiased_moment1 + (moment11 - moment12)
    unbiased_moment2 <- unbiased_moment2 + (moment21 - moment22)
  }
  while(!meeting && (iteration < 1e2)){
    iteration <- iteration + 1
    ## advance two chains with same proposal
    initparticlesprop <- initdist$generate(smctuning$nparticles)
    smc_resultsprop <- smc_hmc_nonadaptive(smctuning, targetdist, initdist$eval, initparticlesprop)
    nstepsprop <- length(smc_resultsprop$lambdas)
    moment1prop <- smc_resultsprop$moment1_history[[nstepsprop]]
    moment2prop  <- smc_resultsprop$moment2_history[[nstepsprop]]
    logzhatprop <- sum(smc_resultsprop$log_ratio_normconst)
    ##
    u_ <- runif(1)
    if (log(u_) < (logzhatprop - logzhat1)){
      ## accept
      meeting <- TRUE
      meetingtime <- iteration
      logzhat1 <- logzhatprop
      ## we can return an unbiased estimator
    } else {
      if (log(u_) < (logzhatprop - logzhat2)){
        logzhat2 <- logzhatprop
        moment12 <- moment1prop
        moment22 <- moment2prop
      }
      ## add bias correction term
      unbiased_moment1 <- unbiased_moment1 + (moment11 - moment12)
      unbiased_moment2 <- unbiased_moment2 + (moment21 - moment22)
    }
  }
  return(list(unbiased_moment1 = unbiased_moment1, unbiased_moment2 = unbiased_moment2,
              meetingtime = meetingtime))
}

#
dim(trainingset)
p <- dim(trainingset)[2]
#
## Prior
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)
priordist <- get_mvnormal_diag(b[,1], diag(B))

##
# load(file = "experiments/logistic/covtype.smcpartial.RData")
# smc_partial_results[[1]]$xmeans_history
# path.partial.df <- data.frame()
# for (rep in 1:length(smc_partial_results)){
#   nsteps <- length(smc_partial_results[[rep]]$xmeans_history)
#   path.partial.df <- rbind(path.partial.df, lapply(1:nsteps, function(time){
#     data.frame(component = 1:p,
#                mean = smc_partial_results[[rep]]$xmeans_history[[time]],
#                var = smc_partial_results[[rep]]$xvars_history[[time]],
#                rep = rep,
#                time = time,
#                ndata = (time - 1) * delta)
#   }) %>% bind_rows())
# }
# tail(path.partial.df)
# ggplot(path.partial.df, aes(x = mean, y = var, group = interaction(component, rep))) +
#   geom_path() + scale_y_log10() + geom_rangeframe()

####
smctuning <- list(nparticles = 2^7, ess_criterion = 0.5, nmoves = 1)
smctuning$stepsize <- 0.3 * p^{-1/4}
smctuning$nleapfrog <- 3
load("experiments/logistic/covtype.smcpartial.RData")
###
ndataseq <- c(0, 25, 50, 100, 250, 500)
ndata_ <- ndataseq[2]

initmean <- b[,1]
initvar <- B

##
ndataseq
nrep <- rep(1e3, 6)
##
unbiased_results <- foreach(irep = 1:nrep[2]) %dorng% {
  subsample <- sample(1:dim(trainingset)[1], size = ndata_, replace = FALSE)
  Y <- trainingset[subsample,1] - 1
  X <- trainingset[subsample,2:11]
  ## standardization of inputs
  X <- apply(X, 2, function(v) (v - mean(v))/(sd(v)))
  X <- cbind(1, X)
  targetdist <- function(xs,...){
    pr <- priordist(xs)
    ll <- smcsamplers:::logistic_loglikelihood_gradient(xs, Y, X)
    return(list(logpdf = pr$logpdf + ll$logls, gradlogpdf = pr$gradlogpdf + ll$gradients))
  }
  initdist <- get_mvnormal(initmean, initvar)
  # initparticles <- initmean + sqrt(initvar) * matrix(rnorm(p*smctuning$nparticles), nrow = p)
  initparticles <- initdist$generate(smctuning$nparticles)
  asmc_results <- asmc_hmc(smctuning, targetdist, initdist$eval, initparticles)
  ##
  variances <- lapply(asmc_results$xhistory, function(x) cov(t(x)))
  smctuning$lambdas <- asmc_results$lambdas
  smctuning$variances <- variances
  result_ <- unbiased_moments(initdist, targetdist, smctuning)
  result_$nsteps <- length(smctuning$lambdas)
  result_
}

sapply(unbiased_results, function(x) x$meetingtime)
sapply(unbiased_results, function(x) x$nsteps)

umean_history <- list()
uvar_history <- list()

# apply(sapply(unbiased_results, function(x) x$unbiased_moment1), 1, sd)
moment1 <- rowMeans(sapply(unbiased_results, function(x) x$unbiased_moment1))
moment2 <- apply(simplify2array(lapply(unbiased_results, function(x) x$unbiased_moment2)), c(1,2), mean)
variance <- moment2 - matrix(moment1, ncol = 1) %*% t(moment1)

initmean <- moment1
initvar <- variance

nboot <- 20

results.df <- data.frame()
for (iboot in 1:nboot){
  results.df <- rbind(results.df, data.frame(component = 1:p,
                                             mean = b[,1],
                                             var =  diag(B),
                                             ndata = 0,
                                             rep = iboot))
}
for (iboot in 1:nboot){
  bootstrapped_ <- unbiased_results[sample(1:nrep[2], size = nrep[2], replace = TRUE)]
  moment1 <- rowMeans(sapply(bootstrapped_, function(x) x$unbiased_moment1))
  moment2 <- apply(simplify2array(lapply(bootstrapped_, function(x) x$unbiased_moment2)), c(1,2), mean)
  variance <- moment2 - matrix(moment1, ncol = 1) %*% t(moment1)
  results.df <- rbind(results.df, data.frame(component = 1:p,
             mean = moment1,
             var =  diag(variance),
             ndata = ndataseq[2],
             rep = iboot))
}
tail(results.df)

for (iobs in 3:length(ndataseq)){
  ndata_ <- ndataseq[iobs]
  cat("# obs:", ndata_, "\n")
  unbiased_results <- foreach(irep = 1:nrep[iobs]) %dorng% {
    subsample <- sample(1:dim(trainingset)[1], size = ndata_, replace = FALSE)
    Y <- trainingset[subsample,1] - 1
    X <- trainingset[subsample,2:11]
    X <- apply(X, 2, function(v) (v - mean(v))/(sd(v)))
    X <- cbind(1, X)
    ##
    targetdist <- function(xs,...){
      pr <- priordist(xs)
      ll <- smcsamplers:::logistic_loglikelihood_gradient(xs, Y, X)
      return(list(logpdf = pr$logpdf + ll$logls, gradlogpdf = pr$gradlogpdf + ll$gradients))
    }
    initdist <- get_mvnormal(initmean, initvar)
    initparticles <- initdist$generate(smctuning$nparticles)
    asmc_results <- asmc_hmc(smctuning, targetdist, initdist$eval, initparticles)
    ##
    variances <- lapply(asmc_results$xhistory, function(x) cov(t(x)))
    smctuning$lambdas <- asmc_results$lambdas
    smctuning$variances <- variances
    result_ <- unbiased_moments(initdist, targetdist, smctuning)
    result_$nsteps <- length(smctuning$lambdas)
    result_
  }
  moment1 <- rowMeans(sapply(unbiased_results, function(x) x$unbiased_moment1))
  moment2 <- apply(simplify2array(lapply(unbiased_results, function(x) x$unbiased_moment2)), c(1,2), mean)
  initmean <- moment1
  initvar <- moment2 - matrix(moment1, ncol = 1) %*% t(moment1)
  # umean_history[[iobs]] <- initmean
  # uvar_history[[iobs]] <- initvar
  print(sapply(unbiased_results, function(x) x$meetingtime))
  print(sapply(unbiased_results, function(x) x$nsteps))
  for (iboot in 1:nboot){
    bootstrapped_ <- unbiased_results[sample(1:nrep[iobs], size = nrep[iobs], replace = TRUE)]
    moment1 <- rowMeans(sapply(bootstrapped_, function(x) x$unbiased_moment1))
    moment2 <- apply(simplify2array(lapply(bootstrapped_, function(x) x$unbiased_moment2)), c(1,2), mean)
    variance <- moment2 - matrix(moment1, ncol = 1) %*% t(moment1)
    results.df <- rbind(results.df, data.frame(component = 1:p,
                                               mean = moment1,
                                               var =  diag(variance),
                                               ndata = ndataseq[iobs],
                                               rep = iboot))
  }
}

# nsteps <- length(umean_history)

# upath.df <- lapply(1:nsteps, function(time){
#   data.frame(component = 1:p,
#              mean = umean_history[[time]],
#              var = diag(uvar_history[[time]]),
#              time = time,
#              ndata = ndataseq[time])
# }) %>% bind_rows()

results.df %>% tail

save(results.df, smctuning, ndataseq, nrep, file = "experiments/logistic/covtype.partial.average.RData")

ggplot(results.df, aes(x = mean, y = var, group = interaction(component, rep))) +
  geom_path() + scale_y_log10() + geom_rangeframe()

# ggplot(path.partial.df, aes(x = mean, y = var, group = interaction(component, rep))) +
#   geom_path() + scale_y_log10() + geom_rangeframe() +
#   geom_path(data=upath.df, aes(x = mean, y = var, group = interaction(component)),
#             colour = 'red')
# load("experiments/logistic/covtype.smctempering.RData")
# path.hmc.df <- data.frame()
# for (rep in 1:length(smc_hmc_results)){
#   nsteps <- length(smc_hmc_results[[rep]]$lambdas)
#   path.hmc.df <- rbind(path.hmc.df, lapply(1:nsteps, function(time){
#     data.frame(component = 1:p,
#                mean = smc_hmc_results[[rep]]$xmeans_history[[time]],
#                var = smc_hmc_results[[rep]]$xvars_history[[time]],
#                rep = rep,
#                time = time)
#   }) %>% bind_rows())
# }
#
# ggplot(path.hmc.df, aes(x = mean, y = var, group = interaction(component, rep))) +
#   geom_path() + scale_y_log10() + geom_rangeframe() +
#   geom_path(data=upath.df, aes(x = mean, y = var, group = interaction(component)),
#             colour = 'red')
#
# load("experiments/logistic/covtype.laplace.n1e+05.RData")
# xlaplace <- smc_results_laplace[[1]]$xhistory
# xlaplace <- xlaplace[[length(xlaplace)]]
# rowMeans(xlaplace)
# diag(cov(t(xlaplace)))
# ggplot(upath.df, aes(x = mean, y = var, group = interaction(component))) +
#   geom_path() + scale_y_log10() + geom_rangeframe() +
#   geom_point(data = data.frame(mean = rowMeans(xlaplace),
#                                var = diag(cov(t(xlaplace)))*100), aes(group = NULL))

