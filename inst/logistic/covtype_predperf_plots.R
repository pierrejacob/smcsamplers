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
