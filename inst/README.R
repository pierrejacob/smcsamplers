## This document describes how to reproduce the figures of the article.

## For the genealogical trees (in Section 4), execute the script trees/treeplots.R
## It requires tidyverse, igraph, ggraph
source("inst/trees/treeplots.R")

## For the results with multivariate Normals in Section 3
source("inst/mvnorm/mvnorm_scalingdimension_run.R")
source("inst/mvnorm/mvnorm_moreintermediatesteps_run.R")
source("inst/mvnorm/mvnorm_plots.R")
## The first two scripts generate RData files, the last one generates the plots.

## (note the script mvnorm/mvnorm_asmc_validity.R is a draft tutorial,
## that could be turned into a Rmarkdown document)

## For the covtype data set,
## first obtain data from
## https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary.html
## https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/covtype.libsvm.binary.scale.bz2
## unzip to obtain file covtype.libsvm.binary.scale
## and execute logistic/covtype_process.R
## you might have to change a path ("~/Downloads/" might not be where your downloaded files go).

source("inst/logistic/covtype_process.R") # creates covtype.processed.RData

## For the "moment plot" of Figure 1
source("inst/logistic/covtype_smc_tempering.R") # creates covtype.smctempering.RData
source("inst/logistic/covtype_smc_partial.R") # creates covtype.smcpartial.RData
source("inst/logistic/covtype_smc_pgg.R") # creates covtype.smcpgg.RData
## and then the plots can be generated with
source("inst/logistic/covtype_pathcomparison_plot.R")

## For the other figures related to the covtype data
source("inst/logistic/covtype_priormerging.R") # creates covtype.merging.RData
source("inst/logistic/covtype_predictiveperf.R") # creates covtype.predperf.RData
source("inst/logistic/covtype_laplace.R") # creates covtype.laplace.n*.RData
### and then to  create plots
source("inst/logistic/covtype_asymptotics_plot.R")

## for the "averaged partial posteriors"
## see script
source("inst/logistic/covtype_averagepartial.R") # creates covtype.partial.average.RData
##  and plots are created with
source("inst/logistic/covtype_averagepartial_plot.R")

## For the SIR model,
## the code requires 'outbreaks', 'rstan'
source("inst/sirmodel/boardingschool-smc-partialposteriors.R") # creates boardingschool-smc-partial.RData
source("inst/sirmodel/boardingschool-marginal.R") # creates boardingschool-smc-marginal.RData
## the plots are created with
source("inst/sirmodel/boardingschool-plots.R")
