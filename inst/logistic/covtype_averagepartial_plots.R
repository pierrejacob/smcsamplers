rm(list = ls())
library(smcsamplers)
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
colours <- c("black", "cornflowerblue", "antiquewhite3")

load(file = "experiments/logistic/covtype.partial.average.RData")
linesize <- 0.05

gpath <- ggplot(results.df, aes(x = mean, y = var, group = interaction(component, rep))) +
  geom_path(alpha = 0.1, size = linesize) + geom_rangeframe()
gpath <- gpath + geom_path(size = linesize, colour = "antiquewhite3") + scale_y_log10()
gpath <- gpath + xlab("mean") + ylab("variance")
gpath <- gpath + scale_x_continuous(breaks = c(-1,0,1))
gpath <- gpath + geom_point(aes(colour = ndata))
gpath <- gpath + scale_color_gradient2(name = "# obs", midpoint = 250,
                                             mid = rgb(.2,0.05,0.2), low = colours[3], high = colours[2], breaks = c(0,500))
gpath <- gpath +  theme(legend.key.width = unit(1.,"cm"))
gpath <- gpath + guides(colour = guide_colorbar(title.vjust=1))
gpath
ggsave(filename = "experiments/logistic/path.partialaverage.pdf", plot = gpath,
       width = 8, height = 8)



