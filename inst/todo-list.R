## Script Figure 1: the Wasserstein distance between two Normals has an explicit formula, why not use it?

## Script Figure 2: it's in matlab!

## smc.R has main function run_smc but in fact refers to a tempered scheme
## and has non-generic feature e.g.
## if(mcmc_kernel$choice == "slice_gibbs")

## ... it might be worth re-thinking about the entire thing

## also it stores the whole history in xtrajectory instead of progressively dumping it to hard-drive
## using e.g. netcdf -> investigate that, would be useful to know anyway and could be a good
## topic for a blog post...

## function get_mvnormal uses some funky terminology: rinit actually samples from the Normal
## so could be called rmvnormal or something like that

## package should not use multiple packages for mvnorm
## either mvtnorm, or mvnfast, or something custom (probably the best option as it would reduce package dependencies)

## figure out netcdf stuff... see test_netcdf.R

## figure out variance estimation... see test_varestim.R

## Cox process example... see test_coxprocess.R

## MVN example... see test_adaptivesmc_hmcmoves.R and test_adaptivesmc_rwmhmoves.R

