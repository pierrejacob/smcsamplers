rm(list = ls())
library(smcsamplers)
library(ggplot2)
set.seed(1)
## gg ridges for the evolution of marginals
theme_set(theme_tufte(ticks = TRUE))
theme_update(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
             axis.title.x = element_text(size = 25, margin = margin(20, 0, 0, 0), hjust = 1),
             axis.title.y = element_text(size = 25, angle = 90, margin = margin(0, 20, 0, 0), vjust = 1),
             legend.text = element_text(size = 20),
             legend.title = element_text(size = 20), title = element_text(size = 30),
             strip.text = element_text(size = 25), strip.background = element_rect(fill = "white"),
             legend.position = "bottom")

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
