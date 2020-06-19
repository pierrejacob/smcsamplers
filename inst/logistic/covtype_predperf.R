rm(list = ls())
library(smcsamplers)
set.seed(1)
library(tidyverse)
library(ggplot2)
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
library(doParallel)
library(doRNG)
registerDoParallel(cores = 6)
load("experiments/logistic/covtype.processed.RData")
#
dim(trainingset)
#
n <- 1e4
ntest <- 1e4
p <- dim(trainingset)[2]-1
Ytest <- testset[1:ntest,1] - 1
Xtest <- testset[1:ntest,2:(p+1)]
Y <- trainingset[1:n,1] - 1
X <- trainingset[1:n,2:(p+1)]

## standardization of inputs
X <- apply(X, 2, function(v) (v - mean(v))/(sd(v)))
X <- cbind(1, X)
Xtest <- apply(Xtest, 2, function(v) (v - mean(v))/(sd(v)))
Xtest <- cbind(1, Xtest)
p <- p + 1

dim(testset)

## adaptive SMC sampler for a chunk of size delta
asmc_hmc_chunk <- function(smctuning, particles, logweights, roots, nstart, delta, Y, X, b, B){
  priordist <- get_mvnormal_diag(b[,1], diag(B))
  initdist <- function(xs,...){
    pr <- priordist(xs)
    if (nstart > 0){
      ll <- smcsamplers:::logistic_loglikelihood_gradient(xs, Y[1:nstart], X[1:nstart,,drop=F])
      return(list(logpdf = pr$logpdf + ll$logls, gradlogpdf = pr$gradlogpdf + ll$gradients))
    } else {
      return(list(logpdf = pr$logpdf, gradlogpdf = pr$gradlogpdf))
    }
  }
  n_next <- min(length(Y), nstart + delta)
  targetdist <- function(xs,...){
    pr <- priordist(xs)
    ll <- smcsamplers:::logistic_loglikelihood_gradient(xs, Y[1:n_next], X[1:n_next,,drop=F])
    return(list(logpdf = pr$logpdf + ll$logls, gradlogpdf = pr$gradlogpdf + ll$gradients))
  }
  particles$init <- initdist(particles$x)
  particles$target <- targetdist(particles$x)
  ## initialize the inverse temperature
  lambda_current <- 0
  lambdas <- c(lambda_current)
  ## initialize the log normalizing constant estimator
  log_ratio_normconst <- c()
  ess_realized <- c()
  infos_mcmc <- list()
  while(lambda_current < 1){
    ## find the next inverse temperature
    logliks <- smcsamplers:::logistic_loglikelihood_gradient(particles$x, Y[(nstart+1):(n_next)], X[(nstart+1):(n_next),])
    ess_deltalambda <- function(lambda){
      logwincremental <- (lambda - lambda_current) * logliks$logls
      logw <- logweights + logwincremental
      return(1/sum(smcsamplers::normalize_weight(logw)$nw^2))
    }
    if (ess_deltalambda(1) > smctuning$ess_criterion * smctuning$nparticles){
      lambda_next <- 1
    } else {
      search_results <- smcsamplers::search_lambda(lambda_current, ess_deltalambda, smctuning$ess_criterion * smctuning$nparticles)
      lambda_next <- search_results$x
    }
    ##
    logwincremental <- (lambda_next - lambda_current) * logliks$logls
    logwincremental[is.na(logwincremental)] <- -Inf
    logweights <- logweights + logwincremental
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
    ## compute effective sample size
    ess_realized <- c(ess_realized, 1/sum(nweights$nw^2))
    ## compute normalizing constant ratio estimator
    log_ratio_normconst <- c(log_ratio_normconst, nweights$avew)
    ##
    lambda_current <- lambda_next
    lambdas <- c(lambdas, lambda_current)
    if (lambda_current < 1){
      ## multinomial resampling
      ancestors <- sample(x = 1:particles$n, size = particles$n,
                          prob = nweights$nw, replace = TRUE)
      ##
      roots <- roots[ancestors]
      particles$x <- particles$x[,ancestors,drop=F]
      particles$init$logpdf <- particles$init$logpdf[ancestors]
      particles$init$gradlogpdf <- particles$init$gradlogpdf[,ancestors,drop=F]
      particles$target$logpdf <- particles$target$logpdf[ancestors]
      particles$target$gradlogpdf <- particles$target$gradlogpdf[,ancestors,drop=F]
      logweights <- rep(0, particles$n)
      ## MCMC move at current lambda
      # accept_rate_ <- 0
      # for (imove in 1:smctuning$nmoves){
      # hmc_res <- hmc_kernel_chunk(particles, smctuning, Y[1:n_next], X[1:n_next,], delta, lambda_current)
      # particles <- hmc_res$particles
      # accept_rate_ <- accept_rate_ + hmc_res$accept_rate
      # }
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
      infos_mcmc[[length(infos_mcmc)+1]] <- info
      # cat("assimilating data from", nstart, "to", nstart+delta, "; current lambda =", lambda_current, "\n")
      # cat("HMC mean acceptance", accept_rate_ / smctuning$nmoves, "\n")
    }
  }
  return(list(particles = particles, logweights = logweights, lambdas = lambdas,
              log_ratio_normconst = log_ratio_normconst,
              ess_realized = ess_realized, infos_mcmc = infos_mcmc, roots = roots))
}

## SMC sampler assimilating data by chunks
## of size delta
asmc_hmc_partial <- function(smctuning, delta, Y, X, b, B){
  starttime <- Sys.time()
  ndata <- length(Y)
  p <- dim(X)[2]
  ## start particles from prior
  priordist <- get_mvnormal_diag(b[,1], diag(B))
  initparticles <- b[,1] + sqrt(diag(B)) * matrix(rnorm(p*smctuning$nparticles), nrow = p)
  particles <- list()
  particles$n <- smctuning$nparticles
  particles$x <- initparticles
  particles$d <- dim(particles$x)[1]
  # particles$init <- initdist(particles$x)
  # particles$target <- targetdist(particles$x)
  logweights <- rep(0, particles$n)
  ##
  istep <- 1
  xmeans_history <- list()
  xvars_history <- list()
  xmeans_history[[1]] <- rowMeans(particles$x)
  xvars_history[[1]] <- apply(particles$x, 1, var)
  nassimilated <- 0
  nassimilatedhistory <- c(nassimilated)
  log_ratio_normconst <- c()
  roots <- 1:particles$n
  nroots <- particles$n
  ##
  xhistory <- list()
  xhistory[[1]] <- particles$x
  whistory <- list()
  whistory[[1]] <- rep(1./particles$n, particles$n)
  lambdas_list <- list()
  infos_mcmc_list <- list()
  while(nassimilated < ndata){
    print(nassimilated)
    asmc_hmc_chunk_results <- asmc_hmc_chunk(smctuning, particles, logweights, roots, nassimilated, delta, Y, X, b, B)
    logweights <- asmc_hmc_chunk_results$logweights
    particles <- asmc_hmc_chunk_results$particles
    roots <- asmc_hmc_chunk_results$roots
    nroots <- c(nroots, length(unique(roots)))
    lambdas_list[[length(lambdas_list)+1]] <- asmc_hmc_chunk_results$lambdas
    infos_mcmc_list[[length(lambdas_list)+1]] <- asmc_hmc_chunk_results$infos_mcmc
    log_ratio_normconst <- c(log_ratio_normconst, sum(asmc_hmc_chunk_results$log_ratio_normconst))
    nassimilated <- min(ndata, nassimilated + delta)
    istep <- istep + 1
    nweights <- smcsamplers::normalize_weight(logweights)
    estimated_means <- sapply(1:particles$d, function(component)
      sum(particles$x[component,] * nweights$nw))
    estimated_variances <- sapply(1:particles$d, function(component)
      sum((particles$x[component,]-estimated_means[component])^2 * nweights$nw))
    xhistory[[istep]] <- particles$x
    whistory[[istep]] <- smcsamplers::normalize_weight(logweights)$nw
    xmeans_history[[istep]] <- estimated_means
    xvars_history[[istep]]  <- estimated_variances
    nassimilatedhistory <- c(nassimilatedhistory, nassimilated)
  }
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(log_ratio_normconst = log_ratio_normconst, nassimilatedhistory = nassimilatedhistory,
              elapsedtime = elapsedtime, xmeans_history = xmeans_history, xvars_history = xvars_history,
              nroots = nroots, lambdas_list = lambdas_list, infos_mcmc_list = infos_mcmc_list,
              xhistory = xhistory, whistory = whistory))
}



smctuning <- list(nparticles = 2^10, ess_criterion = 0.5, nmoves = 2)
smctuning$stepsize <- 0.3 * p^{-1/4}
smctuning$nleapfrog <- ceiling(1/smctuning$stepsize)
delta <- 100
ndataseq <- seq(from = 0, to = length(Y), by = delta)

## Prior
b <- matrix(0, nrow = p, ncol = 1)
B <- diag(10, p, p)

nrep <- 10
smc_partial_results <- foreach(irep = 1:nrep) %dorng% {
  asmc_hmc_partial(smctuning, delta, Y, X, b, B)
}

# save(n, delta, Y, X, b, B, smctuning, smc_partial_results,
#      file = paste0("experiments/logistic/covtype.partial.long.n", n, ".RData"))

## predictive performance
nrep <- length(smc_partial_results)
perf.df <- foreach (rep = 1:nrep, .combine = rbind) %dopar% {
  results <- smc_partial_results[[rep]]
  nsteps <- length(results$xhistory)
  performances_ <- c()
  for (time in 1:nsteps){
    xs_ <- results$xhistory[[time]]
    ws_ <- results$whistory[[time]]
    logliketest <- smcsamplers:::logistic_loglikelihood_gradient(xs_, Ytest, Xtest)$logls
    ## predictive performance:
    ## log integrate p(ytest | beta) posterior(d beta)
    ## log (sum_{n=1}^N W_n p(ytest | beta_n))
    ## log {C * sum_{n=1}^N W_n exp(log p(ytest | beta_n) - log(C))}
    ## log C + log { sum_{n=1}^N W_n exp(log p(ytest | beta_n) - log(C)) }
    # mean(logliketest)
    max_ll_test <- max(logliketest)
    perf <- max_ll_test + log(sum(ws_ * exp(logliketest - max_ll_test)))
    performances_ <- c(performances_, perf)
  }
  data.frame(steps = 1:nsteps,
             nassimilated = (0:(nsteps-1))*delta, performances_, rep = rep)
}
head(perf.df)
save(smctuning, perf.df, ntest, n, file = "experiments/logistic/covtype.predperf.RData")

gperf1 <- ggplot(perf.df %>% filter(nassimilated >= 100, nassimilated <= 5000),
                 aes(x = nassimilated, y = performances_, group = rep))
gperf1 <- gperf1 + geom_line(size = 0.1)
gperf1 <- gperf1 + xlab("# observations") + ylab("log score on test")
gperf1

gperf2 <- ggplot(perf.df %>% filter(nassimilated >= 5000), aes(x = nassimilated, y = performances_, group = rep))
gperf2 <- gperf2 + geom_line(size = 0.1)
gperf2 <- gperf2 + xlab("# observations") + ylab("log score on test")
gperf2
