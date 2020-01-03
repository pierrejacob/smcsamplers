#
set.seed(1)
rm(list = ls())
## test of functions
# temperingsmc:::rmvnorm_
# temperingsmc:::rmvnorm_cholesky_()
# temperingsmc:::dmvnorm_()
# temperingsmc:::dmvnorm_cholesky_inverse()
# temperingsmc:::grad_dmvnorm_precision()

dimension <- 2
rho <- 0.8
# targetdist_mean <- rep(1, dimension)
targetdist_mean <- c(0.2, 0.9)
targetdist_variance <- matrix(0, nrow = dimension, ncol = dimension)
targetdist_variance[1,1] <- 2
targetdist_variance[1,2] <- targetdist_variance[2,1] <- -1
targetdist_variance[2,2] <- 1.5

targetdist_cholvariance <- chol(targetdist_variance)
targetdist_cholinvvariance <- t(chol(solve(targetdist_variance)))

##
# xx <- temperingsmc:::rmvnorm_(nsamples = 1e5, targetdist_mean, targetdist_variance)
# colMeans(xx)
# cov(xx)

##
# xx <- temperingsmc:::rmvnorm_cholesky_(nsamples = 1e5, targetdist_mean, targetdist_cholvariance)
# colMeans(xx)
# cov(xx)

##
xx <- temperingsmc:::rmvnorm_cholesky_(nsamples = 3, targetdist_mean, targetdist_cholvariance)
mvtnorm::dmvnorm(xx, targetdist_mean, targetdist_variance, log = TRUE)
temperingsmc:::dmvnorm_(xx, targetdist_mean, targetdist_variance)
temperingsmc:::dmvnorm_cholesky_inverse(xx, targetdist_mean, targetdist_cholinvvariance)

##
numDeriv::grad(function(x) temperingsmc:::dmvnorm_(x, targetdist_mean, targetdist_variance), xx[1,,drop=F])
numDeriv::grad(function(x) temperingsmc:::dmvnorm_(x, targetdist_mean, targetdist_variance), xx[2,,drop=F])
numDeriv::grad(function(x) temperingsmc:::dmvnorm_(x, targetdist_mean, targetdist_variance), xx[3,,drop=F])
##
temperingsmc:::grad_dmvnorm_precision(xx, targetdist_mean, solve(targetdist_variance))
##
