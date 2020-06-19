## script to test various functions relating to SMC samplers
## initial distribution and target distribution are both multivariate Normals
## geometric bridge between these distributions
## SMC samplers is run with adaptive selection of inverse temperature lambda and resampling
## HMC moves are performed with adaptation of the mass matrix

rm(list = ls())
library(smcsamplers)
library(ggplot2)
graphsettings <- set_custom_theme()
library(doParallel)
library(doRNG)
library(tidyverse)
registerDoParallel(cores = 4)
set.seed(3)
library(ggthemes)
theme_set(theme_tufte(ticks = TRUE))
theme_update(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
             axis.title.x = element_text(size = 25, margin = margin(20, 0, 0, 0), hjust = 1),
             axis.title.y = element_text(size = 25, angle = 90, margin = margin(0, 20, 0, 0), vjust = 1),
             legend.text = element_text(size = 20),
             legend.title = element_text(size = 20), title = element_text(size = 30),
             strip.text = element_text(size = 25), strip.background = element_rect(fill = "white"),
             legend.position = "bottom")

## set up problem
dimension <- 5
initmean <- rep(-0.2, dimension)
initvar <- rep(0.3, dimension)
targetmean <- rep(1.2, dimension)
targetvar <- rep(1.4, dimension)
initdist <- get_mvnormal_diag(initmean, initvar)
targetdist <- get_mvnormal_diag(targetmean, targetvar)

## get exact mean and variance of intermediate distributions
meanvar_intermediate <- function(lambda){
  precision <- (1-lambda) * (1/initvar) + lambda * (1/targetvar)
  var <- 1/precision
  mean <- var * ((1-lambda) * (1/initvar) * initmean + lambda * (1/targetvar) * targetmean)
  return(list(mean = mean, var = var))
}

## set tuning parameters for SMC sampler
smctuning <- list()
## number of particles
smctuning$nparticles <- 2^10
## ESS criterion (number in (0,1))
smctuning$ess_criterion <- 0.75
## stepsize for the HMC moves
smctuning$stepsize <- 1 * dimension^{-1/4}
## number of leapfrog steps
smctuning$nleapfrog <- ceiling(1/smctuning$stepsize)
## number of HMC moves to perform
smctuning$nmoves <- 2
## one run:
initparticles <- initmean + sqrt(initvar) * matrix(rnorm(dimension*smctuning$nparticles), nrow = dimension)
results <- asmc_hmc(smctuning, targetdist, initdist, initparticles)
## what's in the results?
names(results)
nsteps <- length(results$xhistory)

## function to put all x-particles in a data frame
## for the specified components of the state space
get_x <- function(results, components){
  xhistory <- results$xhistory
  lapply(1:length(xhistory), function(time){
    df_ <- data.frame(t(xhistory[[time]])[,components])
    names(df_) <- paste0("X", components)
    df_ %>% mutate(particle = 1:smctuning$nparticles,
                                               lambda = results$lambdas[time],
                                               time = time) %>%
      pivot_longer(cols = starts_with("X"), names_to = "component")
  }
  ) %>% bind_rows()
}

## first, verify validity of sampler
## get estimated means of first component
gmeans <- qplot(results$lambdas, y = sapply(1:nsteps, function(time){
  mean(results$xhistory[[time]][1,])
}), geom = 'line')
gmeans <- gmeans + xlab(expression(lambda)) + ylab(expression(E(X[1])))
gmeans <- gmeans + geom_line(aes(y = sapply(results$lambdas, function(lambda) meanvar_intermediate(lambda)$mean[1])),
                             colour = 'red')
##
## get estimated variances of first component
gvars <- qplot(results$lambdas, y = sapply(1:nsteps, function(time){
  var(results$xhistory[[time]][1,])
}), geom = 'line')
gvars <- gvars + xlab(expression(lambda)) + ylab(expression(V(X[1])))
gvars <- gvars + geom_line(aes(y = sapply(results$lambdas, function(lambda) meanvar_intermediate(lambda)$var[1])),
                           colour = 'red')
gridExtra::grid.arrange(gmeans, gvars, nrow = 2)

## plot two marginals at some intermediate step
## overlaid with exact marginal pdf functions
df_ <- get_x(results, c(3,5))
ilambda <- floor(nsteps/2)
lambda <- results$lambdas[ilambda]
g <- ggplot(df_ %>% filter(time == ilambda), aes(x = value)) + geom_histogram(aes(y=..density..))
g <- g + facet_wrap(~component)
mv_ <- meanvar_intermediate(lambda)
g <- g + stat_function(fun = function(x) dnorm(x, mean = mv_$mean[1], sd = sqrt(mv_$var[1])), colour = 'red')
g <- g + labs(subtitle = paste("\\lambda =", round(lambda, 2)))
g <- g + geom_rangeframe()
print(g)

## ridge plot for the evolution of the marginal
## distributions of a component (e.g. the first one)
library(ggridges)
comp <- (df_$component %>% unique)[1]
gridges <- ggplot() + ylim(0,1.1)
# g <- g + geom_density_ridges(alpha = 0.25, stat = "density", position = "identity")
gridges <- gridges + xlab(expression(x[1])) + ylab(expression(lambda))
ridgescale <- 0.25
for (ilambda in 1:length(results$lambdas)){
  lambda_ <- results$lambdas[ilambda]
  mv <- meanvar_intermediate(lambda_)
  m_ <- mv$mean[1]
  s_ <- sqrt(mv$var[1])
  df_ridge <- data.frame(x = seq(from = m_-3*s_, to = m_+3*s_, length.out = 100))
  df_ridge$y = lambda_
  df_ridge$height = dnorm(df_ridge$x, mean = m_, sd = s_)
  x1s <- results$xhistory[[ilambda]][1,]
  nw <- results$nwhistory[[ilambda]]
  f_ <- density(x = x1s, weights = nw)
  df_hist <- data.frame(x = f_$x, height = f_$y, y = lambda_)
  gridges <- gridges + geom_ridgeline(data=df_hist, aes(x = x, y = y, height = height, group = NULL), alpha = 0,
                          col = 'red', fill = 'red', scale = ridgescale)
  gridges <- gridges + geom_ridgeline(data=df_ridge, aes(x = x, y = y, height = height, group = NULL), alpha = 0.,
                          scale = ridgescale)
}
gridges <- gridges + coord_flip()
print(gridges)

## get information about how the MCMC steps performed
info_mcmc_df <- lapply(2:length(results$lambdas), function(ilambda) {
  info_mcmc <- results$infos_mcmc[[ilambda-1]]
  df_ <- data.frame(ilambda = ilambda, imove = 1:smctuning$nmoves,
                    ar = sapply(info_mcmc, function(x) x$ar),
                    sqjd = sapply(info_mcmc, function(x) x$sqjd))
}) %>% bind_rows()
##
gar <- ggplot(info_mcmc_df, aes(x = ilambda, y = ar)) + geom_point()
gar <- gar + xlab("step") + ylab("acceptance rate") + ylim(0,1) + geom_rangeframe()
gsqjd <- ggplot(info_mcmc_df, aes(x = ilambda, y = sqjd)) + geom_point()
gsqjd <- gsqjd + xlab("step") + ylab("relative jumping distance") + geom_rangeframe()
gridExtra::grid.arrange(gar, gsqjd, nrow = 2)
##

## plot genealogical tree
genealogy_ <- ahistory2genealogy(results$ahistory)
library(igraph)
library(ggraph)
head(genealogy_$dendro)
## uncomment following to create graphs of genealogy
mygraph <- graph_from_data_frame(genealogy_$dendro)
ggraph(mygraph, layout = 'dendrogram', circular = F) +
  geom_edge_diagonal() +
  theme_void() + coord_flip() + scale_y_reverse()

## plot of nunber of unique ancestors
guniqueroots <- qplot(x = results$lambdas, y = sapply(results$roots_history, function(x) length(unique(x))))
guniqueroots <- guniqueroots + xlab(expression(lambda)) + ylab("# unique roots") + ylim(0, results$particles$n)
guniqueroots

### next create some plots and calculations related to multiple independent runs
## multiple runs
smctuning$nparticles <- 2^9
smctuning$nmoves <- 1
nrep <- 100
indep_results <- foreach(irep = 1:nrep) %dorng% {
  initparticles <- initmean + sqrt(initvar) * matrix(rnorm(dimension*smctuning$nparticles), nrow = dimension)
  asmc_hmc(smctuning, targetdist, initdist, initparticles)
}

## plot marginal of some components at final time
df_ <- lapply(1:50, function(r){
  get_x(indep_results[[r]], c(3,5)) %>% mutate(rep = r) %>% filter(time == max(time))
}) %>% bind_rows()
g <- ggplot(df_, aes(x = value, group = rep)) + geom_density(aes(y=..density..))
g <- g + facet_wrap(~component)
mv_ <- meanvar_intermediate(1)
g <- g + stat_function(fun = function(x) dnorm(x, mean = mv_$mean[1], sd = sqrt(mv_$var[1])), colour = 'red')
g <- g + geom_rangeframe()
print(g)

## number of intermediate steps
summary(sapply(indep_results, function(x) length(x$xhistory)))

## Z hat
zhat <- sapply(indep_results, function(x) sum(x$log_ratio_normconst))
summary(zhat)
var(zhat)

## variance of Zhat
var_estims <- sapply(indep_results, function(x) variance_estimator(x))
meanZ <- mean(exp(zhat))
aggreg_var_estim <- mean((exp(zhat)/meanZ)^2 * var_estims)
aggreg_var_estim
var(zhat)

## plot number of roots as time progresses
uniquerootsdf <- lapply(1:length(indep_results), function(x){
  data.frame(rep = x,
             lambdas = indep_results[[x]]$lambdas,
             nuniqueroots = sapply(indep_results[[x]]$roots_history, function(v) length(unique(v))))
  }) %>% bind_rows()
guniqueroots <- ggplot(uniquerootsdf, aes(x = lambdas, y = nuniqueroots, group = rep)) + geom_point()
guniqueroots <- guniqueroots + xlab(expression(lambda)) + ylab("# roots") + ylim(0, smctuning$nparticles)
guniqueroots + geom_rangeframe()


