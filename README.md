# aSPU2 

aSPU2 is a new version for aSPU and related tests. We rewrite the codes in a more efficient way and expect aSPU2 is around 10 times faster than the previous verison (aSPU). Please cite the following manuscript for primary aSPU method:

* Pan, Wei, et al. "A powerful and adaptive association test for rare variants." Genetics 197.4 (2014): 1081-1095.

## Features

* Very efficient codes to calculate the p-values. It is around ten times faster than the original aSPU package.
* It's more general than the aSPU package. The input is score vector and its corresponding covariance matrix.
* Incorporate external weights, improving statistical power.
* Two new ways to combine multiple weights (aSPUO and daSPU).


## Installation

To install the stable version from CRAN, run the following from an R console: 

```R
install.packages("aSPU2")
```

To install the latest development builds directly from GitHub, run this instead:

```R
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("ChongWu-Biostat/aSPU2")
```

For some servers, we have to specify the package directory and install the package locally. "/gpfs/home/cwu3/R/x86_64-redhat-linux-gnu-library/3.5/" is the directory that stores the R package, and you need to replace it with your own directory.
```R
library(devtools)
library(withr)
withr::with_libpaths(new = "/gpfs/home/cwu3/R/x86_64-redhat-linux-gnu-library/3.5/", install_github("ChongWu-Biostat/aSPU2"))
```

*Please make sure OpenMP is installed successfully. To speed up, the core function in aSPU is written in C++ and some software is needed. See [link](http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/#after-3-4-0) for details. *



### How to use

* First, we generate score, covariance matrix and some weights.

```R
library(aSPU2)
n = 100
times <- 1:n
rho <- 0.3
sigma <- 2

###############
H <- abs(outer(times, times, "-"))
V <- sigma * rho^H
p <- nrow(V)
V[cbind(1:p, 1:p)] <- V[cbind(1:p, 1:p)] * sigma

eV<-eigen(V)
CovSsqrt<- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))

U = CovSsqrt %*% rnorm(n)

# Multiple weights we want to use
weight = cbind(rep(1,100),runif(100))
```

* Second, we run aSPU2.

```R
# with single weight
aSPU(U,V,weight[,1],pow=c(1:6,Inf),n.perm = 10000)

```



## Stay In Touch

- For latest releases and announcements, follow on GitHub: [@ChongWu-Biostat](https://github.com/ChongWu-Biostat)


## License

[MIT](http://opensource.org/licenses/MIT)

Copyright (c) 2013-present, Chong Wu (cwu3@fsu.edu) 

