rm(list = ls())
library(smcsamplers)
library(ggplot2)
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

colours <- c("black", "cornflowerblue", "antiquewhite3")

load("experiments/logistic/covtype.smctempering.RData")
linesize <- 0.05
plotwidth = 5
plotheight = 8

##
p <- smc_hmc_results[[1]]$particles$d
path.hmc.df <- data.frame()
for (rep in 1:length(smc_hmc_results)){
  nsteps <- length(smc_hmc_results[[rep]]$lambdas)
  path.hmc.df <- rbind(path.hmc.df, lapply(1:nsteps, function(time){
    data.frame(component = 1:p,
    mean = smc_hmc_results[[rep]]$xmeans_history[[time]],
    var = smc_hmc_results[[rep]]$xvars_history[[time]],
    rep = rep,
    time = time,
    lambda = smc_hmc_results[[rep]]$lambdas[time])
  }) %>% bind_rows())
}
path.hmc.df[path.hmc.df$time==1,2] <- 0
path.hmc.df[path.hmc.df$time==1,3] <- 10
gpathtempering <- ggplot(path.hmc.df, aes(x = mean, y = var, group = interaction(component, rep)))
gpathtempering <- gpathtempering + geom_path(size = linesize, colour = "antiquewhite3") + scale_y_log10()
gpathtempering <- gpathtempering + xlab("mean") + ylab("variance")
gpathtempering <- gpathtempering + scale_x_continuous(breaks = c(-1,0,1))
gpathtempering <- gpathtempering + geom_point(aes(colour = lambda))
# gpathtempering <- gpathtempering + viridis::scale_color_viridis(name = expression(lambda), discrete=F, breaks = c(0,1)) #+ theme(legend.position = 'none')
gpathtempering <- gpathtempering + scale_color_gradient2(name = expression(lambda), midpoint = 0.5,
                     mid = rgb(.2,0.05,0.2), low = colours[3], high = colours[2], breaks = c(0,1))
gpathtempering <- gpathtempering +  theme(legend.key.width = unit(1.,"cm"))
gpathtempering <- gpathtempering + guides(colour = guide_colorbar(title.vjust=1))
gpathtempering

ggsave(filename = "experiments/logistic/path.tempering.pdf", plot = gpathtempering, width = plotwidth, height = plotheight)

load("experiments/logistic/covtype.smcpgg.RData")
##
p <- asmc_pgg_results[[1]]$particles$d
path.pgg.df <- data.frame()
for (rep in 1:length(asmc_pgg_results)){
  nsteps <- length(asmc_pgg_results[[rep]]$lambdas)
  path.pgg.df <- rbind(path.pgg.df, lapply(1:nsteps, function(time){
    data.frame(component = 1:p,
               mean = asmc_pgg_results[[rep]]$xmeans_history[[time]],
               var = asmc_pgg_results[[rep]]$xvars_history[[time]],
               rep = rep,
               time = time,
               lambda = asmc_pgg_results[[rep]]$lambdas[time])
  }) %>% bind_rows())
}

# ggplot(path.pgg.df, aes(x = mean, y = var, group = interaction(component, rep))) + geom_path() + scale_y_log10()
path.pgg.df[path.pgg.df$time==1,2] <- 0
path.pgg.df[path.pgg.df$time==1,3] <- 10

gpathpgg <- ggplot(path.pgg.df, aes(x = mean, y = var, group = interaction(component, rep)))
gpathpgg <- gpathpgg + geom_path(size = linesize, colour = "antiquewhite3") + scale_y_log10()
gpathpgg <- gpathpgg + xlab("mean") + ylab("variance")
gpathpgg <- gpathpgg + scale_x_continuous(breaks = c(-1,0,1))
gpathpgg <- gpathpgg + geom_point(aes(colour = lambda))
gpathpgg <- gpathpgg + scale_color_gradient2(name = expression(lambda), midpoint = 0.5,
                           mid = rgb(.2,0.05,0.2), low = colours[3], high = colours[2], breaks = c(0,1))
# +  viridis::scale_color_viridis(name = expression(lambda), discrete=F, breaks = c(0,1)) #+ theme(legend.position = 'none')
gpathpgg <- gpathpgg +  theme(legend.key.width = unit(1.,"cm"))
gpathpgg <- gpathpgg + guides(colour = guide_colorbar(title.vjust=1))
gpathpgg
ggsave(filename = "experiments/logistic/path.pgg.pdf", plot = gpathpgg, width = plotwidth, height = plotheight)

load("experiments/logistic/covtype.smcpartial.RData")

path.partial.df <- data.frame()
for (rep in 1:length(smc_partial_results)){
  nsteps <- length(smc_partial_results[[rep]]$xmeans_history)
  path.partial.df <- rbind(path.partial.df, lapply(1:nsteps, function(time){
    data.frame(component = 1:p,
               mean = smc_partial_results[[rep]]$xmeans_history[[time]],
               var = smc_partial_results[[rep]]$xvars_history[[time]],
               rep = rep,
               time = time,
               ndata = smc_partial_results[[rep]]$nassimilatedhistory[time])
  }) %>% bind_rows())
}
path.partial.df[path.partial.df$time==1,2] <- 0
path.partial.df[path.partial.df$time==1,3] <- 10

ggplot(path.partial.df, aes(x = mean, y = var, group = interaction(component, rep))) +
  geom_path() + scale_y_log10() + geom_rangeframe()
gpathpartial <- ggplot(path.partial.df, aes(x = mean, y = var, group = interaction(component, rep)))
gpathpartial <- gpathpartial + geom_path(size = linesize, colour = "antiquewhite3") + scale_y_log10()
gpathpartial <- gpathpartial + xlab("mean") + ylab("variance")
gpathpartial <- gpathpartial + scale_x_continuous(breaks = c(-1,0,1))
gpathpartial <- gpathpartial + geom_point(aes(colour = ndata))

gpathpartial <- gpathpartial + scale_color_gradient2(name = "# obs", midpoint = 500,
                                             mid = rgb(.2,0.05,0.2), low = colours[3], high = colours[2], breaks = c(0,1000))
# +  viridis::scale_color_viridis(name = "# obs", discrete=F, breaks = c(0,1000)) #+ theme(legend.position = 'none')
gpathpartial <- gpathpartial +  theme(legend.key.width = unit(1.,"cm"))
gpathpartial <- gpathpartial + guides(colour = guide_colorbar(title.vjust=1))
gpathpartial
ggsave(filename = "experiments/logistic/path.partial.pdf", plot = gpathpartial, width = plotwidth, height = plotheight)



# ## agreement on final time
# ggplot(path.pgg.df %>% filter(time == max(time)),
#        aes(x = component, y = mean)) + geom_point(position = position_jitter()) +
#   geom_point(data = path.hmc.df %>% filter(time == max(time)), colour = 'red', alpha = 0.5, position = position_jitter()) +
#   geom_point(data = path.partial.df %>% filter(time == max(time)), colour = 'blue', alpha = 0.5, position = position_jitter())
#
# ggplot(path.pgg.df %>% filter(time == max(time)),
#        aes(x = component, y = var)) + geom_point(position = position_jitter()) +
#   geom_point(data = path.hmc.df %>% filter(time == max(time)), colour = 'red', alpha = 0.5, position = position_jitter()) +
#   geom_point(data = path.partial.df %>% filter(time == max(time)), colour = 'blue', alpha = 0.5, position = position_jitter()) +
#   scale_y_log10()

# ## performance of MCMC moves
# ## get information about how the MCMC steps performed
# hmcperf.df <- data.frame()
# for (rep in 1:length(smc_hmc_results)){
#   info_mcmc_df <- lapply(2:length(smc_hmc_results[[rep]]$lambdas), function(ilambda) {
#     info_mcmc <- smc_hmc_results[[rep]]$infos_mcmc[[ilambda-1]]
#     df_ <- data.frame(ilambda = ilambda, imove = 1:length(info_mcmc),
#                       sqjd = sapply(info_mcmc, function(x) x$sqjd))
#   }) %>% bind_rows()
#   hmcperf.df <- rbind(hmcperf.df, info_mcmc_df %>% mutate(rep = rep))
# }
# ##
# gsqjd <- ggplot(hmcperf.df, aes(x = ilambda, y = sqjd)) + geom_point()
# gsqjd <- gsqjd + xlab("step") + ylab("relative jumping distance") + geom_rangeframe()
# gsqjd
#
# pggperf.df <- data.frame()
# for (rep in 1:length(smc_hmc_results)){
#   info_mcmc_df <- lapply(2:length(asmc_pgg_results[[rep]]$lambdas), function(ilambda) {
#     info_mcmc <- asmc_pgg_results[[rep]]$infos_mcmc[[ilambda-1]]
#     df_ <- data.frame(ilambda = ilambda, imove = 1:length(info_mcmc),
#                       sqjd = sapply(info_mcmc, function(x) x$sqjd))
#   }) %>% bind_rows()
#   pggperf.df <- rbind(pggperf.df, info_mcmc_df %>% mutate(rep = rep))
# }
# ##
# gsqjd <- ggplot(pggperf.df, aes(x = ilambda, y = sqjd)) + geom_point()
# gsqjd <- gsqjd + xlab("step") + ylab("relative jumping distance") + geom_rangeframe()
# gsqjd


