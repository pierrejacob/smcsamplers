
rm(list = ls())
library(smcsamplers)
library(ggplot2)
library(ggthemes)
graphsettings <- set_custom_theme()
library(doParallel)
library(doRNG)
library(tidyverse)
registerDoParallel(cores = 6)
set.seed(3)
theme_set(theme_tufte(ticks = TRUE))
theme_update(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
             axis.title.x = element_text(size = 25, margin = margin(20, 0, 0, 0), hjust = 1),
             axis.title.y = element_text(size = 25, angle = 90, margin = margin(0, 20, 0, 0), vjust = 1),
             legend.text = element_text(size = 20),
             legend.title = element_text(size = 20), title = element_text(size = 30),
             strip.text = element_text(size = 25), strip.background = element_rect(fill = "white"),
             legend.position = "bottom")

load(file = "experiments/mvnorm/scalingfixedn.RData")
head(fixedndf)

load(file = "experiments/mvnorm/scalinglinearn.RData")
head(linearndf)

load(file = "experiments/mvnorm/scalingintermediate.RData")
head(intermedf)

colours <- c("black", "cornflowerblue", "antiquewhite3")

##
# gnsteps <- ggplot(fixedndf, aes(x = dimension, y = nsteps, group = rep)) + geom_point()
# gsqerror <- ggplot(fixedndf  %>% group_by(dimension) %>% summarize(meansqerror = mean(sqerror)),
#                    aes(x = dimension, y = meansqerror)) + geom_point()
# gzhat <- ggplot(fixedndf %>% group_by(dimension) %>% summarise(varzhat = var(zhat)),
#                 aes(x = dimension, y = varzhat)) + geom_point()
# gnroots <- ggplot(fixedndf, aes(x = dimension, y = nroots, group = rep)) + geom_point()
# gridExtra::grid.arrange(gnsteps, gnroots, gsqerror, gzhat)
#
# gnsteps <- ggplot(linearndf, aes(x = dimension, y = nsteps, group = rep)) + geom_point()
# gsqerror <- ggplot(linearndf  %>% group_by(dimension) %>% summarize(meansqerror = mean(sqerror)),
#                    aes(x = dimension, y = meansqerror)) + geom_point()
# gzhat <- ggplot(linearndf %>% group_by(dimension) %>% summarise(varzhat = var(zhat)),
#                 aes(x = dimension, y = varzhat)) + geom_point()
# gnroots <- ggplot(linearndf, aes(x = dimension, y = nroots, group = rep)) + geom_point()
# gridExtra::grid.arrange(gnsteps, gnroots, gsqerror, gzhat)
#
# gnsteps <- ggplot(intermedf, aes(x = dimension, y = nsteps, group = rep)) + geom_point()
# gsqerror <- ggplot(intermedf  %>% group_by(dimension) %>% summarize(meansqerror = mean(sqerror)),
#                    aes(x = dimension, y = meansqerror)) + geom_point()
# gzhat <- ggplot(intermedf %>% group_by(dimension) %>% summarise(varzhat = var(zhat)),
#                 aes(x = dimension, y = varzhat)) + geom_point()
# gnroots <- ggplot(intermedf, aes(x = dimension, y = nroots, group = rep)) + geom_point()
# gridExtra::grid.arrange(gnsteps, gnroots, gsqerror, gzhat)


df_ <- rbind(fixedndf %>% mutate(case = "fixed N"),
  linearndf %>% mutate(case = "linear N"),
  intermedf %>% mutate(case = "fixed N & d steps"))
head(df_)

plotwidth <- 7
plotheight <- 5
linesize <- 1
gnsteps <- ggplot(df_ %>% filter(case != "fixed N & d steps") %>% group_by(dimension, case) %>% summarise(meannsteps = mean(nsteps)),
       aes(x = dimension, y = meannsteps, group = interaction(case), colour = case, linetype = case))
gnsteps <- gnsteps + geom_line(size = linesize)
gnsteps <- gnsteps + scale_color_manual(name = '', values = colours[1:2]) + ylab("# steps")
gnsteps <- gnsteps + scale_linetype_discrete(name = "") + geom_rangeframe()
gnsteps <- gnsteps + scale_x_continuous(breaks = dimensions)
gnsteps
ggsave(filename = "experiments/mvnorm/nsteps.pdf", plot = gnsteps, width = plotwidth, height = plotheight)

gnroots <- ggplot(df_ %>% group_by(dimension, case) %>% summarise(meannroots = mean(nroots)),
                  aes(x = dimension, y = meannroots, colour = case, linetype = case)) + geom_line(size = linesize)
gnroots <- gnroots + scale_color_manual(name = '', values = colours[c(1,3,2)]) + ylab("mean # roots")
gnroots <- gnroots + scale_linetype(name = "")
gnroots <- gnroots + geom_rangeframe() + scale_x_continuous(breaks = dimensions)
gnroots
ggsave(filename = "experiments/mvnorm/nroots.pdf", plot = gnroots, width = plotwidth, height = plotheight)

gsqerror <- ggplot(df_  %>% group_by(dimension, case) %>% summarize(meansqerror = mean(sqerror)),
                   aes(x = dimension, y = meansqerror, colour = case, linetype = case)) + geom_line(size = linesize)
gsqerror <- gsqerror + scale_color_manual(name = '', values = colours[c(1,3,2)]) + ylab("MSE E[X]") + scale_linetype(name = "")
gsqerror <- gsqerror + geom_rangeframe()  + scale_x_continuous(breaks = dimensions)
gsqerror
ggsave(filename = "experiments/mvnorm/sqerror.pdf", plot = gsqerror, width = plotwidth, height = plotheight)

gzhat <- ggplot(df_ %>% group_by(dimension, case) %>% summarise(varzhat = var(zhat)),
                aes(x = dimension, y = varzhat, colour = case, linetype = case)) + geom_line(size = linesize)
gzhat <- gzhat + scale_color_manual(name = '', values = colours[c(1,3,2)]) + ylab("variance log constant")  + scale_linetype(name = "")
gzhat <- gzhat + geom_rangeframe() + scale_x_continuous(breaks = dimensions)
gzhat
ggsave(filename = "experiments/mvnorm/varzhat.pdf", plot = gzhat, width = plotwidth, height = plotheight)

