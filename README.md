
## smcsamplers

Implements SMC samplers to re-create the figures of the article entitled
*An invitation to sequential Monte Carlo samplers*, written by Chenguang
Dai, Jeremy Heng, Pierre E. Jacob, and Nick Whiteley, [available on
arXiv at arxiv.org/abs/2007.11936](https://arxiv.org/abs/2007.11936).

This is not a general-purpose statistical software. This is just a
collection of scripts intended to reproduce the figures of a paper. Use
at your own risk!

The folder inst/README.R describes which scripts generate which figures.

### Installation

The package can be installed from R via:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("pierrejacob/smcsamplers")
```

It depends on the packages `Rcpp`, `RcppEigen`, `tidyverse`,
`doParallel`, `doRNG`, `ggplot2`, `ggridges`, `ggthemes`, which are all
on CRAN.
