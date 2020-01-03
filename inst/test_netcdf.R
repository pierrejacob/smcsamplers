# install.packages('ncdf4')
# library(ncdf4)
rm(list = ls())
set.seed(1)

## file name
ncfname <- "~/test_ncdf4.nc"
## create two matrices, a vector and a real
nsamples <- 23
dimension <- 5
samples1 <- matrix(rnorm(nsamples * dimension), nrow = nsamples, ncol = dimension)
samples2 <- matrix(rnorm(nsamples * dimension), nrow = nsamples, ncol = dimension)
nyy <- 17
yy <- rgamma(nyy, 3, 1)
xx <- rgamma(1, 3, 3)
##
## we will create a netcdf file with samples1
## and then add samples2, yy and xx one by one
##
## create dimensions
ncdim_samples1nrow <- ncdf4::ncdim_def("samples1nrow", "count", vals = 1:nsamples)
ncdim_samples1ncol <- ncdf4::ncdim_def("samples1ncol", "count", vals = 1:dimension)
## create variable
ncvar_samples1 <- ncdf4::ncvar_def(name = "samples1", units = "blurp", dim = list(ncdim_samples1nrow, ncdim_samples1ncol))

## create netcdf file with samples1 in it
nc_new <- ncdf4::nc_create(filename = ncfname, vars = ncvar_samples1)
ncdf4::ncvar_put(nc = nc_new, varid = "samples1", vals = samples1)
ncdf4::nc_close(nc_new)

## add vector yy to file
ncdim_yy <- ncdf4::ncdim_def("yylength", "count", vals = 1:nyy)
ncvar_yy <- ncdf4::ncvar_def(name = 'yy', units = 'blurp', dim = ncdim_yy)

nc_connect <- ncdf4::nc_open(filename = ncfname, write = TRUE)
ncdf4::ncvar_add(nc = nc_connect, v = ncvar_yy)
ncdf4::nc_close(nc_connect)

nc_connect <- ncdf4::nc_open(filename = ncfname, write = TRUE)
ncdf4::ncvar_put(nc = nc_connect, 'yy', yy)
ncdf4::nc_close(nc_connect)

## add samples2 to file
## samples2 has same size as samples1 so no need to create new dimensions
ncvar_samples2 <- ncdf4::ncvar_def(name = "samples2", units = "blurp", dim = list(ncdim_samples1nrow, ncdim_samples1ncol))
nc_connect <- ncdf4::nc_open(filename = ncfname, write = TRUE)
ncdf4::ncvar_add(nc = nc_connect, v = ncvar_samples2)
ncdf4::nc_close(nc_connect)

nc_connect <- ncdf4::nc_open(filename = ncfname, write = TRUE)
ncdf4::ncvar_put(nc = nc_connect, 'samples2', samples2)
ncdf4::nc_close(nc_connect)

## finally add a scalar
ncvar_xx <- ncdf4::ncvar_def(name = "xx", units = "blurp", dim = list())
nc_connect <- ncdf4::nc_open(filename = ncfname, write = TRUE)
ncdf4::ncvar_add(nc = nc_connect, v = ncvar_xx)
ncdf4::nc_close(nc_connect)
## add value
nc_connect <- ncdf4::nc_open(filename = ncfname, write = TRUE)
ncdf4::ncvar_put(nc_connect, varid = 'xx', vals = xx)
ncdf4::nc_close(nc_connect)


## Now get these variables back
nc_connect <- ncdf4::nc_open(filename = ncfname)
for (ivar in 1:length(nc_connect$var)){
  varname = nc_connect$var[[ivar]]$name
  print(varname)
}

retrieved_yy <- ncdf4::ncvar_get(nc_connect, varid = "yy")
retrieved_samples1 <- ncdf4::ncvar_get(nc_connect, varid = "samples1")
retrieved_samples2 <- ncdf4::ncvar_get(nc_connect, varid = "samples2")
retrieved_xx <- ncdf4::ncvar_get(nc_connect, varid = "xx")
ncdf4::nc_close(nc_connect)


summary(retrieved_yy - yy)
summary(retrieved_samples1 - samples1)
summary(retrieved_samples2 - samples2)
summary(retrieved_xx - xx)

### add large matrices recursively to a netcdf file

# niterations <- 10
# # history_ <- list()
# # for (iteration in 1:niterations){
# #   history_[[iteration]] <- rnorm(5e6)
# # }
# nsamples <- 5e6
# ncfname <- "~/test_large_ncdf4.nc"
# ncdim_x <- ncdf4::ncdim_def("xlength", "count", vals = 1:nsamples)
# ncvar_ <- ncdf4::ncvar_def(name = paste0('x', 1), units = 'blurp', dim = ncdim_x)
# nc_connect <- ncdf4::nc_create(filename = ncfname, vars = ncvar_)
# ncdf4::ncvar_put(nc = nc_connect, varid = "x1", vals = rnorm(nsamples))
# ncdf4::nc_close(nc_connect)
#
# for (iteration in 2:niterations){
#   nc_connect <- ncdf4::nc_open(filename = ncfname, write = TRUE)
#   ncvar_ <- ncdf4::ncvar_def(name = paste0('x', iteration), units = 'blurp', dim = ncdim_x)
#   ncdf4::ncvar_add(nc = nc_connect, v = ncvar_)
#   ncdf4::nc_close(nc_connect)
#   nc_connect <- ncdf4::nc_open(filename = ncfname, write = TRUE)
#   ncdf4::ncvar_put(nc_connect, varid = paste0('x', iteration), vals = rnorm(nsamples))
#   ncdf4::nc_close(nc_connect)
# }

