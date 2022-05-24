rm(list = ls())
library(smcsamplers)
graphsettings <- set_custom_theme()
theme_set(theme_tufte(ticks = TRUE))
theme_update(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
             axis.title.x = element_text(size = 25, margin = margin(20, 0, 0, 0), hjust = 1),
             axis.title.y = element_text(size = 25, angle = 90, margin = margin(0, 20, 0, 0), vjust = 1),
             legend.text = element_text(size = 20),
             legend.title = element_text(size = 20), title = element_text(size = 30),
             strip.text = element_text(size = 25), strip.background = element_rect(fill = "white"),
             legend.position = "bottom")
set.seed(3) # for reproductibility
library(outbreaks)
library(tidyverse)
colours <- c("black", "cornflowerblue", "antiquewhite3")

# time series of cases
cases <- influenza_england_1978_school$in_bed # Number of students in bed
flu <- data.frame(day = seq_along(cases), cases = cases)
gflu <- ggplot(flu, aes(x = day, y = cases)) + geom_point()
gflu <- gflu + xlab("days") + ylab("# confined in bed") #+ labs(subtitle = "# confined in bed", title = "Influenza in a boarding school")
gflu <- gflu + scale_x_continuous(breaks = c(3,7,11,14))
gflu <- gflu + scale_y_continuous(breaks = 0:3 * 100)
gflu <- gflu + ggthemes::geom_rangeframe()
gflu

ggsave(filename = "experiments/sirmodel/boardingschool.flu.pdf", plot = gflu, width = 7, height = 5)

load("experiments/sirmodel/boardingschool-smc-partial.RData")

xhistory <- smc_partial_results[[4]]
length(xhistory)

df_hist <- data.frame()
# for (time in 1:2){
for (ndata in 1:length(xhistory)){
  for (component in 1:3){
    x1s <- xhistory[[ndata]][component,]
    f_ <- density(x = x1s)
    df_hist <- rbind(df_hist, data.frame(x = f_$x, height = f_$y, ndata = ndata + 2,
                                         component = component))
  }
}
df_hist %>% tail

ridgescale <- 1
gridges <- ggplot()
gridges <- gridges + xlab(expression(phi[inv])) + ylab("# observations")
gridges <- gridges + geom_ridgeline(data=df_hist %>% filter(component == 3), aes(x = exp(x), y = ndata, height = height, group = ndata),
                                    alpha = 0.5, col = colours[2], fill = colours[2], size = 0.2)
gridges <- gridges + scale_y_continuous(breaks = c(3,7,11,14)) + xlim(0,1)
gridges <- gridges + coord_flip()
gridges
ggsave(filename = "experiments/sirmodel/boardingschool.partial.phiinv.pdf", plot = gridges, width = 7, height = 5)

##
df_moments <- lapply(1:length(smc_partial_results), function(index){
  xhistory <- smc_partial_results[[index]]
  df_moments_ <- data.frame()
  for (ndata in 1:length(xhistory)){
    xs <- xhistory[[ndata]]
    df_moments_ <- rbind(df_moments_, data.frame(component = 1:3,
                                               means = rowMeans(xs),
                                               vars = apply(xs, 1, var),
                                               ndata = ndata+2,
                                               rep = index))
  }
  df_moments_
}) %>% bind_rows()

head(df_moments)

gmoments <- ggplot(df_moments,aes(x = means, y = vars, group = interaction(component, rep))) + geom_path()
gmoments <- gmoments + scale_y_log10()
gmoments

xhistory <- smc_partial_results[[1]]
xs <- xhistory[[1]]
dim(xs)
plot(xs[1,], xs[2,], col = rainbow(14)[1])
for (ndata in 2:length(xhistory)){
  xs <- xhistory[[ndata]]
  points(xs[1,], xs[2,], col = rainbow(14)[ndata])
}

pointsdf <- lapply(1:length(xhistory), function(ndata){
  data.frame(t(xhistory[[ndata]])) %>% mutate(ndata = ndata + 2)
}) %>% bind_rows()

gpoints <- ggplot(pointsdf, aes(x = exp(X1), y = exp(X2), colour = ndata)) + geom_point() +
  scale_color_gradient2(name = "# observations", midpoint = mean(pointsdf$ndata),
                        low = "black", mid = "white", high = "cornflowerblue",
                        breaks = c(3,7,11,14))
gpoints <- gpoints + xlab(expression(gamma)) + ylab(expression(beta))
gpoints <- gpoints + theme(legend.key.width = unit(1.5,"cm"))
gpoints <- gpoints + guides(col = guide_colorbar(nbin = 5))
gpoints
# ggsave(filename = "experiments/sirmodel/boardingschool.partial.points.pdf", plot = gpoints, width = 7, height = 5)

gcont <- ggplot(pointsdf %>% filter(ndata == 3), aes(x = exp(X1), y = exp(X2), group = ndata)) +
  stat_density_2d(aes(fill = ndata), bins = 2, geom = 'polygon', colour = 'black', size = 0.1) +
  stat_density_2d(data = pointsdf %>% filter(ndata == 4), aes(fill = ndata), bins = 2, geom = 'polygon', colour = 'black', size = 0.1) +
  stat_density_2d(data = pointsdf %>% filter(ndata == 5), aes(fill = ndata), bins = 2, geom = 'polygon', colour = 'black', size = 0.1) +
  stat_density_2d(data = pointsdf %>% filter(ndata == 6), aes(fill = ndata), bins = 2, geom = 'polygon', colour = 'black', size = 0.1) +
  stat_density_2d(data = pointsdf %>% filter(ndata == 7), aes(fill = ndata), bins = 2, geom = 'polygon', colour = 'black', size = 0.1) +
  stat_density_2d(data = pointsdf %>% filter(ndata == 8), aes(fill = ndata), bins = 2, geom = 'polygon', colour = 'black', size = 0.1) +
  stat_density_2d(data = pointsdf %>% filter(ndata == 9), aes(fill = ndata), bins = 2, geom = 'polygon', colour = 'black', size = 0.1) +
  stat_density_2d(data = pointsdf %>% filter(ndata == 10), aes(fill = ndata), bins = 2, geom = 'polygon', colour = 'black', size = 0.1) +
  stat_density_2d(data = pointsdf %>% filter(ndata == 11), aes(fill = ndata), bins = 2, geom = 'polygon', colour = 'black', size = 0.1) +
  stat_density_2d(data = pointsdf %>% filter(ndata == 12), aes(fill = ndata), bins = 2, geom = 'polygon', colour = 'black', size = 0.1) +
  stat_density_2d(data = pointsdf %>% filter(ndata == 13), aes(fill = ndata), bins = 2, geom = 'polygon', colour = 'black', size = 0.1) +
  stat_density_2d(data = pointsdf %>% filter(ndata == 14), aes(fill = ndata), bins = 2, geom = 'polygon', colour = 'black', size = 0.1)
## I am aware that the above code looks offensive, but couldn't make anything else work
gcont <- gcont +  scale_fill_gradient2(name = "# observations", midpoint = mean(pointsdf$ndata),
                        low = "antiquewhite3", mid = "white", high = "cornflowerblue",
                        breaks = c(3,7,11,14))
gcont <- gcont + xlab(expression(gamma)) + ylab(expression(beta))
gcont <- gcont + theme(legend.key.width = unit(1.5, "cm"))
gcont <- gcont + guides(col = guide_colorbar(nbin = 5))
gcont <- gcont  +  scale_fill_gradient2(name = "# observations", midpoint = mean(pointsdf$ndata),
                       low = "antiquewhite3", mid = "white", high = "cornflowerblue",
                       breaks = c(3,7,11,14))
gcont

ggsave(filename = "experiments/sirmodel/boardingschool.partial.contours.pdf", plot = gcont, width = 7, height = 5)


load(file = "experiments/sirmodel/boardingschool-smc-marginal.RData")
head(resultsdf)
gzhat <- ggplot(resultsdf, aes(x = t0, y = zhat, group = irep)) + geom_line(size = 0.25) + xlab(expression(t[0])) + ylab("Z hat")
gzhat <- gzhat + ylab("log marginal likelihood")
gzhat
ggsave(filename = "experiments/sirmodel/boardingschool.partial.logz.pdf", plot = gzhat, width = 7, height = 5)

