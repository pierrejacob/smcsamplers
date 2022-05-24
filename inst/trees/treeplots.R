rm(list=ls())
set.seed(1)
library(smcsamplers)
if (!require("igraph")) install.packages("igraph")
library(igraph)
if (!require("ggraph")) install.packages("ggraph")
library(ggraph)
if (!require("gridExtra")) install.packages("gridExtra")
graphsettings <- set_custom_theme()
set.seed(3)

## initial distribution is bivariate Normal
dimension <- 2
initmean <- rep(5, dimension)
initvar <- rep(5, dimension)
initdist <- get_mvnormal_diag(initmean, initvar)
## target distribution is "banana"
# target log-density evaluation
bananatarget <- function(x) -(1-x[1])^2 - 10*((x[2]-x[1]^2)^2)
# gradient of target log-density
bananagradtarget <- function(x) c(-2*(x[1]-1) + 40*x[1]*(x[2]-x[1]^2),
                                  -20 * (x[2]-x[1]^2))
targetdist <- function(xs){
  return(list(gradlogpdf = apply(xs, 2, bananagradtarget), logpdf = apply(xs, 2, bananatarget)))
}
## SMC samplers tuning parameters
smctuning <- list()
smctuning$ess_criterion <- 0.5
smctuning$nmoves <- 2
asmc_results <- list()
smctuning$nparticles <- 2^8
smctuning$stepsize <- 0.1
smctuning$nleapfrog <- ceiling(1/smctuning$stepsize)
initparticles <- initmean + sqrt(initvar) * matrix(rnorm(dimension*smctuning$nparticles), nrow = dimension)
## run adaptive SMC sampler
results <- asmc_hmc(smctuning, targetdist, initdist, initparticles)
nsteps <- length(results$lambdas)
print(nsteps)

## get information about how the MCMC steps performed
info_mcmc_df <- lapply(2:length(results$lambdas), function(ilambda) {
  info_mcmc <- results$infos_mcmc[[ilambda-1]]
  df_ <- data.frame(ilambda = ilambda, imove = 1:smctuning$nmoves,
                    ar = sapply(info_mcmc, function(x) x$ar),
                    sqjd = sapply(info_mcmc, function(x) x$sqjd))
}) %>% bind_rows()
##
gar <- ggplot(info_mcmc_df, aes(x = ilambda, y = ar)) + geom_point()
gar <- gar + xlab("step") + ylab("acceptance rate") + ylim(0,1)
gsqjd <- ggplot(info_mcmc_df, aes(x = ilambda, y = sqjd)) + geom_point()
gsqjd <- gsqjd + xlab("step") + ylab("relative jumping distance")
gridExtra::grid.arrange(gar, gsqjd, nrow = 2)

## plot genealogical tree
genealogy_ <- ahistory2genealogy(results$ahistory)
head(genealogy_$dendro)
mygraph <- graph_from_data_frame(genealogy_$dendro)

genealogygraph <- ggraph(mygraph, layout = 'dendrogram', circular = T) +
  geom_edge_diagonal(aes(colour = factor(istep==1)), edge_width = 0.2) +
  theme_void() + scale_edge_color_manual(values = c("black", "white")) +
  theme(legend.position = "none")
## add number of unique ancestors in the middle
genealogygraph <- genealogygraph + annotate(geom = 'label', x = 0, y = 0, label = paste0(length(unique(results$roots_history[[length(results$roots_history)]]))))
print(genealogygraph)
ggsave(filename = "experiments/bigtree.pdf", plot = genealogygraph, width = 6, height = 6)

## plot genealogical tree for SMCS run with smaller number of particles
for (rep in 1:3){
  smctuning$nparticles <- 2^6
  initparticles <- initmean + sqrt(initvar) * matrix(rnorm(dimension*smctuning$nparticles), nrow = dimension)
  results <- asmc_hmc(smctuning, targetdist, initdist, initparticles)
  genealogy_ <- ahistory2genealogy(results$ahistory)
  mygraph <- graph_from_data_frame(genealogy_$dendro)
  g <- ggraph(mygraph, layout = 'dendrogram', circular = F) +
    geom_edge_diagonal(edge_width = 0.2) +
    theme_void() + coord_flip() + scale_y_reverse()
  print(g)
  ggsave(filename = paste0("experiments/smalltree", rep, ".pdf"), plot = g, width = 5, height = 5)
}
