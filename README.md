tcensReg: Maximum Likelihood Estimation of a Truncated Normal Distribution with Censored Data
================

The goal of this package is to estimate parameters from a linear model when the data comes from a truncated normal distribution with censoring. Maximum likelihood values are returned derived from Newton-Raphson algorithm using analytic values of the gradient and hessian. This package is also able to return maximum likelihood estimates for truncated only or censored only data similar to `truncreg` and `censReg` packages.

Example
-------

Some common examples where this type of problem may arise is when there is a natural truncation imposed by the structure of the data. For instance several applications have an implied zero truncation such as product lifetimes, age, or detection thresholds. To show how to implement the functions within the package, I will demonstrate a simple simulation example.

Assume that we have observations from an underlying truncated normal distribution

*Y*<sup>\*</sup> ∼ TN(*μ*, *σ*<sup>2</sup>, *ν*),

where *ν* denotes the value of the left-truncation. In our case we will assume a zero-truncated model by setting *ν* = 0.

``` r
library(msm) #we will use this package to generate random values from the truncated normal distribution
mu <- 0.8
sigma <- 0.5
nu <- 0

y_star <- msm::rtnorm(n = 100, mean = mu, sd = sigma, lower = nu)
range(y_star) #note that the lowerbound will always be non-negative
```

    ## [1] 0.1082646 2.4518047

Next, we can imagine a scenario where we have an imprecise measurement of *Y*<sup>\*</sup> leading to values below *a* being censored. We will set *a* = 0.5 in this example creating the random variable *Y*, where

*Y*...

``` r
a <- 0.5
y <- ifelse(y_star<=a, a, y_star)
sum(y==a)/length(y) #calculating the number of censored observations
```

    ## [1] 0.28

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

``` r
#installing the package for estimating truncated with censoring from the GitHub page
devtools::install_github("williazo/tcensReg") 
```

    ## Skipping install of 'tcensReg' from a github remote, the SHA1 (8248f1dc) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
library(tcensReg)  #loading the package into the current environment
```
