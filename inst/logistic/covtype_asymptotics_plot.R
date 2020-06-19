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

load(file = "experiments/logistic/covtype.merging.RData")
col1 <- "cornflowerblue"
# col1 <- "darkgrey"
col2 <- "antiquewhite3"
ridgescale <- 10
gridges <- ggplot()
gridges <- gridges + xlab(expression(beta[1])) + ylab("# observations") + geom_vline(xintercept = 0, linetype = 2)
gridges <- gridges + geom_ridgeline(data=df_hist %>% filter(component == 1),
                                    aes(x = x, y = ndata, height = y, group = ndata, color = "prior 1", fill = "prior 1"),
                                    scale = ridgescale, alpha  = 0.1, size  = 0.2)
gridges <- gridges + geom_ridgeline(data=df_hist2 %>% filter(component == 1),
                                    aes(x = x, y = ndata, height = y, group = ndata, color = "prior 2", fill = "prior 2"),
                                    scale = ridgescale, alpha  = 0.1, size  = 0.2)
gridges <- gridges + scale_y_continuous(breaks = c(0,50,100,150,200)) + xlim(-5,5)
gridges <- gridges + scale_fill_manual(name = "", labels = c("prior 1", "prior 2"), values = c(col1, col2)) +
  scale_color_manual(name = "", labels = c("prior 1", "prior 2"), values = c(col1, col2))
# gridges <- gridges + guides(color = guid)
print(gridges)
ggsave(filename = "experiments/logistic/covtype.merging.beta1.pdf", plot = gridges, width = 7, height = 5)

ridgescale <- 50
gridges <- ggplot()
gridges <- gridges + xlab(expression(beta[8])) + ylab("# observations") + geom_vline(xintercept = 0, linetype = 2)
gridges <- gridges + geom_ridgeline(data=df_hist %>% filter(component == 8),
                                    aes(x = x, y = ndata, height = y, group = ndata, color = "prior 1", fill = "prior 1"),
                                    scale = ridgescale, alpha  = 0.1, size  = 0.2)
gridges <- gridges + geom_ridgeline(data=df_hist2 %>% filter(component == 8),
                                    aes(x = x, y = ndata, height = y, group = ndata, color = "prior 2", fill = "prior 2"),
                                    scale = ridgescale, alpha  = 0.1, size  = 0.2)
gridges <- gridges + scale_y_continuous(breaks = c(0,50,100,150,200)) + xlim(-5,5)
gridges <- gridges + scale_fill_manual(name = "", labels = c("prior 1", "prior 2"), values = c(col1, col2)) +
  scale_color_manual(name = "", labels = c("prior 1", "prior 2"), values = c(col1, col2))
# gridges <- gridges + guides(color = guid)
print(gridges)

ggsave(filename = "experiments/logistic/covtype.merging.beta8.pdf", plot = gridges, width = 7, height = 5)

# load("experiments/logistic/covtype.partial.long.n1000.RData")
# p <- ncol(X)
#
# results <- smc_partial_results[[1]]
# ## compute ESS against observations
# nsteps <- length(results$xhistory)
# plot((0:(nsteps-1))*delta, sapply(1:nsteps, function(time) 1/sum(results$whistory[[time]]^2)),
#      type = 'l')
# abline(v = ((0:(nsteps-1))*delta)[sapply(1:nsteps, function(time) length(results$infos_mcmc_list[[time]]))>0])
#
# results$infos_mcmc_list
#
# results$nassimilatedhistory
# results$elapsedtime
# results$nroots
# length(results$xhistory)
# length(results$whistory)
# length(results$xmeans_history)
# ## check that xmeans_history can be retrieved from xhistory, whistory
# # time <- 12
# # results$xmeans_history[[time]]
# # sapply(1:p, function(comp) sum(results$whistory[[time]] * results$xhistory[[time]][comp,]))
#
# path.partial.df <- data.frame()
# for (rep in 1:length(smc_partial_results)){
#   nsteps <- length(smc_partial_results[[rep]]$xmeans_history)
#   path.partial.df <- rbind(path.partial.df, lapply(1:nsteps, function(time){
#     data.frame(component = 1:p,
#                mean = smc_partial_results[[rep]]$xmeans_history[[time]],
#                var = smc_partial_results[[rep]]$xvars_history[[time]],
#                rep = rep,
#                time = time,
#                nobservations  = (time-1)*delta)
#   }) %>% bind_rows())
# }
# path.partial.df %>% tail
# ggplot(path.partial.df %>% filter(nobservations >= 0), aes(x = mean, y = var, group = interaction(component, rep))) +
#   geom_path() + scale_y_log10()
#
#
# load("experiments/logistic/covtype.partial.long2.n1000.RData")
#
# ## ridge plot
# ## during the first steps...
# component <- 8
# library(ggridges)
# results <- smc_partial_results[[1]]
# results2 <- smc_partial_results2[[1]]
# gridges <- ggplot()
# gridges <- gridges + xlab(expression(beta[8])) + ylab("# observations") + geom_vline(xintercept = 0, linetype = 2)
# ridgescale <- 50
# for (time in 1:15){
#   x1s <- results$xhistory[[time]][component,]
#   nw <- results$whistory[[time]]
#   x1s2 <- results2$xhistory[[time]][component,]
#   nw2 <- results2$whistory[[time]]
#   f_ <- density(x = x1s, weights = nw)
#   f_2 <- density(x = x1s2, weights = nw2)
#   df_hist <- data.frame(x = f_$x, height = f_$y, y = (time-1) * delta)
#   gridges <- gridges + geom_ridgeline(data=df_hist, aes(x = x, y = y, height = height, group = NULL), alpha = 0.3,
#                                       col = 'grey', fill = 'grey', scale = ridgescale, size = 0.2)
#   df_hist2 <- data.frame(x = f_2$x, height = f_2$y, y = (time-1) * delta)
#   gridges <- gridges + geom_ridgeline(data=df_hist2, aes(x = x, y = y, height = height, group = NULL), alpha = 0.3,
#                                       col = 'orange', fill = 'orange', scale = ridgescale, size = 0.2)
# }
# gridges <- gridges + scale_y_continuous(breaks = delta * (0:(14))) + xlim(-5,5)
# gridges <- gridges + coord_flip()
# gridges
# ggsave(filename = "experiments/logistic/covtype.ridge.2priors.pdf", plot = gridges,
#        width = 10, height = 5)

# ## later on, with more observations ...
# component <- 8
# library(ggridges)
# results <- smc_partial_results[[1]]
# results2 <- smc_partial_results2[[1]]
# gridges <- ggplot()
# gridges <- gridges + xlab(expression(x)) + ylab("# observations") + geom_vline(xintercept = 0, linetype = 2)
# ridgescale <- 50
# for (time in 20:length(results$xhistory)){
#   x1s <- results$xhistory[[time]][component,]
#   nw <- results$whistory[[time]]
#   x1s2 <- results2$xhistory[[time]][component,]
#   nw2 <- results2$whistory[[time]]
#   f_ <- density(x = x1s, weights = nw)
#   f_2 <- density(x = x1s2, weights = nw2)
#   df_hist <- data.frame(x = f_$x, height = f_$y, y = (time-1) * delta)
#   gridges <- gridges + geom_ridgeline(data=df_hist, aes(x = x, y = y, height = height, group = NULL), alpha = 0.3,
#                                       col = 'grey', fill = 'grey', scale = ridgescale, size = 0.2)
#   df_hist2 <- data.frame(x = f_2$x, height = f_2$y, y = (time-1) * delta)
#   gridges <- gridges + geom_ridgeline(data=df_hist2, aes(x = x, y = y, height = height, group = NULL), alpha = 0.3,
#                                       col = 'orange', fill = 'orange', scale = ridgescale, size = 0.2)
# }
# gridges <- gridges + coord_flip()
# print(gridges)
# load("experiments/logistic/covtype.partial.long.n10000.RData")
# load("experiments/logistic/covtype.partial.long2.n10000.RData")
#
# path.partial.df <- data.frame()
# for (rep in 1:length(smc_partial_results)){
#   nsteps <- length(smc_partial_results[[rep]]$xmeans_history)
#   path.partial.df <- rbind(path.partial.df, lapply(1:nsteps, function(time){
#     data.frame(component = 1:p,
#                mean = smc_partial_results[[rep]]$xmeans_history[[time]],
#                var = smc_partial_results[[rep]]$xvars_history[[time]],
#                rep = rep,
#                time = time,
#                nobservations  = (time-1)*delta)
#   }) %>% bind_rows())
# }
#
# path.partial.df2 <- data.frame()
# for (rep in 1:length(smc_partial_results2)){
#   nsteps <- length(smc_partial_results2[[rep]]$xmeans_history)
#   path.partial.df2 <- rbind(path.partial.df2, lapply(1:nsteps, function(time){
#     data.frame(component = 1:p,
#                mean = smc_partial_results2[[rep]]$xmeans_history[[time]],
#                var = smc_partial_results2[[rep]]$xvars_history[[time]],
#                rep = rep,
#                time = time,
#                nobservations  = (time-1)*delta)
#   }) %>% bind_rows())
# }
# path.partial.df2 %>% tail
# gpathpartial <- ggplot(path.partial.df, aes(x = mean, y = var, group = interaction(component, rep))) +
#   geom_path(size = 0.1) +
#   geom_path(size = 0.1, data=path.partial.df2, colour = 'orange') + scale_y_log10()
# gpathpartial <- gpathpartial + xlab("mean") + ylab("variance")
# gpathpartial <- gpathpartial + scale_x_continuous(breaks = c(0,1))
# gpathpartial
# ggsave(filename = "experiments/logistic/path.merging.pdf",
#        plot = gpathpartial, width = 5, height = 5)

# ggplot(path.partial.df %>% filter(nobservations >= 5000), aes(x = mean, y = var, group = interaction(component, rep))) +
#   geom_path() +
#   geom_path(data=path.partial.df2 %>% filter(nobservations >= 5000), colour = 'red') + scale_y_log10()


# ## predictive performance
# nrep <- length(smc_partial_results)
# perf.df <- foreach (rep = 1:nrep, .combine = rbind) %dopar% {
#   results <- smc_partial_results[[rep]]
#   nsteps <- length(results$xhistory)
#   performances_ <- c()
#   for (time in 1:nsteps){
#     xs_ <- results$xhistory[[time]]
#     ws_ <- results$whistory[[time]]
#     logliketest <- smcsamplers:::logistic_loglikelihood_gradient(xs_, Ytest, Xtest)$logls
#     ## predictive performance:
#     ## log integrate p(ytest | beta) posterior(d beta)
#     ## log (sum_{n=1}^N W_n p(ytest | beta_n))
#     ## log {C * sum_{n=1}^N W_n exp(log p(ytest | beta_n) - log(C))}
#     ## log C + log { sum_{n=1}^N W_n exp(log p(ytest | beta_n) - log(C)) }
#     # mean(logliketest)
#     max_ll_test <- max(logliketest)
#     perf <- max_ll_test + log(sum(ws_ * exp(logliketest - max_ll_test)))
#     performances_ <- c(performances_, perf)
#   }
#   data.frame(steps = 1:nsteps,
#              nassimilated = (0:(nsteps-1))*delta, performances_, rep = rep)
# }
# head(perf.df)

load(file = "experiments/logistic/covtype.predperf.RData")
gperf1 <- ggplot(perf.df %>% filter(nassimilated >= 100, nassimilated <= 5000),
                 aes(x = nassimilated, y = performances_, group = rep))
gperf1 <- gperf1 + geom_line(size = 0.1)
gperf1 <- gperf1 + xlab("# observations") + ylab("log score on test")
gperf1
ggsave(filename = "experiments/logistic/covtype.pred1.pdf",
       plot = gperf1, width = 7, height = 5)

gperf2 <- ggplot(perf.df %>% filter(nassimilated >= 5000), aes(x = nassimilated, y = performances_, group = rep))
gperf2 <- gperf2 + geom_line(size = 0.1)
gperf2 <- gperf2 + xlab("# observations") + ylab("log score on test")
gperf2
ggsave(filename = "experiments/logistic/covtype.pred2.pdf",
       plot = gperf2, width = 7, height = 5)

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

# gelapsed <- ggplot(laplace_df, aes(x = n, y = elapsedtime, group = rep)) + geom_line() +
#   ylab("elapsed (seconds)") + xlab("# observations") + scale_x_continuous(breaks = c(1e4, 5e4, 1e5))
# gelapsed <- gelapsed + geom_rangeframe()
# gelapsed
# ggsave(filename = "experiments/logistic/covtype.laplace.elapsed.pdf", plot = gelapsed,
#        width = 7, height = 5)

