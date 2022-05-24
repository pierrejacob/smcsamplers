### logistic regression with covtype data set
rm(list = ls())
library(smcsamplers)
set.seed(1)
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

##
## the following function performs some precomputation useful
## when implementing the Polya-Gamma Gibbs sampler
## 'b' and 'B' refer respectively to the mean and covariance matrix of the prior on 'beta'
## (typically b = (0,...,0) while B = diag(sigma2, ..., sigma2)...)

pgg_precomputation <- function(Y, X, b, B){
  invB <- solve(B)
  invBtimesb <- invB %*% b
  Ykappa <- matrix(Y - rep(0.5, length(Y)), ncol=1)
  XTkappa <- t(X) %*% Ykappa
  KTkappaplusinvBtimesb <- XTkappa + invBtimesb
  return(list(n=nrow(X), p=ncol(X), X=X, Y=Y, b=b, B=B,
              invB=invB, invBtimesb=invBtimesb, KTkappaplusinvBtimesb=KTkappaplusinvBtimesb))
}
## function to obtain mean and variance of 'beta' given all the Polya Gamma variables
## in a vector omega, and given precomputed quantities obtained e.g. via 'pgg_precomputation'
pgg_m_and_sigma <- function(omega, precomputed){
  return(smcsamplers:::pgg_m_sigma_(omega, precomputed$X, precomputed$invB, precomputed$KTkappaplusinvBtimesb))
}

## function that takes a vector 'beta' and precomputed quantities (obtained e.g. via 'pgg_precomputation')
## and performs one step of Polya Gamma Gibbs sampler
pgg_kernel <- function(beta, precomputed){
  zs <- abs(precomputed$X %*% beta)
  w <- pgdraw::pgdraw(1, zs)
  res <- pgg_m_and_sigma(w, precomputed)
  beta <- smcsamplers:::rmvnorm_cholesky_(1, res$m, res$Cholesky)
  return(t(beta))
}


##
# pandl <- get_prior_and_ll(Y, X, B[1,1])
# logprior <- pandl$logprior
rinit <- function(n){
  return(b[,1] + sqrt(diag(B)) * matrix(rnorm(p*n), nrow = p))
}


## path of distribution: from prior to posterior via multiplying X by lambda
## adaptive SMC sampler
asmc_pgg <- function(smctuning){
  starttime <- Sys.time()
  particles <- list()
  particles$n <- smctuning$nparticles
  particles$x <- rinit(particles$n)
  particles$d <- dim(particles$x)[1]
  ### initialize the inverse temperature
  lambda_current <- 0
  lambdas <- c(lambda_current)
  ### initialize the log normalizing constant estimator
  log_ratio_normconst <- c()
  ess_ <- c()
  ## history of particles
  istep <- 1
  # xhistory <- list()
  # xhistory[[istep]] <- particles$x
  roots <- 1:particles$n
  nroots <- c(particles$n)
  xmeans_history <- list()
  xvars_history <- list()
  xmeans_history[[1]] <- rowMeans(particles$x)
  xvars_history[[1]] <- apply(particles$x, 1, var)
  infos_mcmc <- list()
  # while inv temperature is not one...
  while(lambda_current < 1){
    print(lambda_current)
    ## find the next inverse temperature
    xbeta <- X %*% particles$x
    Yxbeta_sums <- colSums(Y * xbeta)
    ess_deltalambda <- function(lambda){
      logw <- (lambda - lambda_current) * Yxbeta_sums - colSums(log(1 + exp(lambda * xbeta)) - log(1 + exp(lambda_current * xbeta)))
      return(1/sum(smcsamplers::normalize_weight(logw)$nw^2))
    }
    if (ess_deltalambda(1) > smctuning$ess_criterion * smctuning$nparticles){
      lambda_next <- 1
    } else {
      search_results <- smcsamplers::search_lambda(lambda_current, ess_deltalambda, smctuning$ess_criterion * smctuning$nparticles)
      lambda_next <- search_results$x
    }
    ## computed incremental log-weights
    logweights <- (lambda_next - lambda_current) * Yxbeta_sums - colSums(log(1 + exp(lambda_next * xbeta)) - log(1 + exp(lambda_current * xbeta)))
    logweights[is.na(logweights)] <- -Inf
    nweights <- smcsamplers::normalize_weight(logweights)
    ess_ <- c(ess_, 1/sum(nweights$nw^2))
    log_ratio_normconst <- c(log_ratio_normconst, nweights$avew)
    estimated_means <- sapply(1:particles$d, function(component)
      sum(particles$x[component,] * nweights$nw))
    estimated_variances <- sapply(1:particles$d, function(component)
      sum((particles$x[component,]-estimated_means[component])^2 * nweights$nw))
    ## multinomial resampling
    ancestors <- sample(x = 1:smctuning$nparticles, size = smctuning$nparticles,
                        prob = nweights$nw, replace = TRUE)
    ## keep track of roots of the genealogical tree
    roots <- roots[ancestors]
    particles$x <- particles$x[,ancestors,drop=F]
    lambda_current <- lambda_next
    lambdas <- c(lambdas, lambda_current)
    ## MCMC move at current lambda
    pgg_precomputed <- pgg_precomputation(Y, lambda_current * X, b, B)
    info <- list()
    for (imove in 1:smctuning$nmoves){
      newparticlesx <- apply(particles$x, 2, function(v) pgg_kernel(v, pgg_precomputed))
      sqjd <- mean(colSums((newparticlesx - particles$x)^2/estimated_variances)/particles$d)
      particles$x <- newparticlesx
      info[[imove]] <- list(sqjd = sqjd)
    }
    infos_mcmc[[istep]] <- info
    istep <- istep + 1
    ## store mean and variance
    xmeans_history[[istep]] <- estimated_means
    xvars_history[[istep]] <- estimated_variances
    nroots <- c(nroots, length(unique(roots)))
  }
  currenttime <- Sys.time()
  elapsedtime <- as.numeric(lubridate::as.duration(lubridate::ymd_hms(currenttime) - lubridate::ymd_hms(starttime)), "seconds")
  return(list(particles = particles, lambdas = lambdas, log_ratio_normconst = log_ratio_normconst,
              ess_ = ess_, elapsedtime = elapsedtime, roots = roots, nroots = nroots,
              xmeans_history = xmeans_history, xvars_history = xvars_history, infos_mcmc = infos_mcmc))
}

##
smctuning <- list(nparticles = 2^9, ess_criterion = 0.5, nmoves = 1)

nrep <- 10
asmc_pgg_results <- foreach(irep = 1:nrep) %dorng% {
  asmc_pgg(smctuning)
}
##
save(smctuning, nrep, asmc_pgg_results, file = "experiments/logistic/covtype.smcpgg.RData")



