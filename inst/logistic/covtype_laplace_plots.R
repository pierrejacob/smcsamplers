rm(list = ls())
library(smcsamplers)
library(ggplot2)
plotwidth <- 7
plotheight <- 5
set.seed(1)
library(dplyr)
## gg ridges for the evolution of marginals
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

library(reshape2)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
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


## next Laplace approximation
laplace.files <- list.files(path = "experiments/logistic/", pattern = "laplace")
laplace.files <- laplace.files[grepl(pattern = "RData", x = laplace.files)]
laplace_df <- lapply(1:length(laplace.files), function(ifile){
  load(file = paste0("experiments/logistic/", laplace.files[ifile]))
  nrep <- length(smc_results_laplace)
  nroots <- sapply(smc_results_laplace, function(x){
    nroots_ <- sapply(x$roots_history, function(z) length(unique(z)))
    nroots_[length(nroots_)]})
  elapsedtime <- sapply(smc_results_laplace, function(x) x$elapsedtime)
  data.frame(rep = 1:nrep,
             ess = sapply(smc_results_laplace, function(x) x$ess_realized),
             logZhat = sapply(smc_results_laplace, function(x) x$log_ratio_normconst),
             nroots = nroots,
             n = n,
             elapsedtime = elapsedtime)
}) %>% bind_rows()

load(file = paste0("experiments/logistic/", laplace.files[1]))
laplace_df %>% head

gess <- ggplot(laplace_df, aes(x = n, y = 100* ess/smctuning$nparticles, group = rep)) + geom_point() +
  ylab("ESS %") + xlab("# observations") + scale_x_continuous(breaks = c(1e4, 5e4, 1e5)) +
  ylim(99, 100) + geom_rangeframe()
gess
ggsave(filename = "experiments/logistic/covtype.laplace.ess.pdf", plot = gess, width = 7, height = 5)

glogzhat <- ggplot(laplace_df, aes(x = n, y = logZhat / n, group = rep)) + geom_line() +
  ylab("log constant / # observations") + xlab("# observations") + scale_x_continuous(breaks = c(1e4, 5e4, 1e5))
glogzhat <- glogzhat + geom_rangeframe()
glogzhat
ggsave(filename = "experiments/logistic/covtype.laplace.logzhat.pdf", plot = glogzhat,
       width = 7, height = 5)

