tcensReg: Maximum Likelihood Estimation of a Truncated Normal
Distribution with Censored Data
================

The goal of this package is to estimate parameters from a linear model
when the data comes from a truncated normal distribution with censoring.
Maximum likelihood values are returned derived from Newton-Raphson
algorithm using analytic values of the gradient and hessian. This
package is also able to return maximum likelihood estimates for
truncated only or censored only data similar to `truncreg` and `censReg`
packages.

# Installation

You can install ggplot.spaghetti from github via the devtools package
with:

``` r
install.packages("devtools")
devtools::install_github("williazo/tcensReg")
#to install the package with the accompaining vignette use the command below
devtools::install_github("williazo/tcensReg", build_opts=c("--no-resave-data", "--no-manual"),
                         build_vignettes=TRUE)
```

# Example 1: Single Population

Some common examples where this type of problem may arise is when there
is a natural truncation imposed by the structure of the data. For
instance several applications have an implied zero truncation such as
product lifetimes, age, or detection thresholds. To show how to
implement the functions within the package, I will demonstrate a simple
simulation example.

Assume that we have observations from an underlying truncated normal
distribution. In our case we will assume a zero-truncated model by
setting a=0. We generate this truncated normal data below and refer to
it as
`ystar`.

``` r
library(msm) #we will use this package to generate random values from the truncated normal distribution
mu <- 0.5
sigma <- 0.5
a <- 0

y_star <- msm::rtnorm(n=1000, mean=mu, sd=sigma, lower=a)
range(y_star) #note that the lowerbound will always be non-negative
```

    ## [1] 0.0002510439 2.2261729774

Next, we can imagine a scenario where we have an imprecise measurement
of `ystar` leading to censoring. In our case we assume that values below
a limit of detection, `nu`, are censored. This creates a random variable
`y`.

In the example below we set our limit of detection as `nu`=0.25.

``` r
nu <- 0.25
y <- ifelse(y_star <= nu, nu, y_star)
sum(y == nu)/length(y) #calculating the number of censored observations
```

    ## [1] 0.187

``` r
dt <- data.frame(y_star, y) #collecting the uncensored and censored data together
```

We can observe the histogram and density plot for the uncensored data,
which shows the zero-truncation.
![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

We can then compare this to the censored observations below
![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

We can then estimate the mean, `mu`, and standard deviation, `sigma`,
using `y` with the `tcensReg` package as shown below.

``` r
library(tcensReg)  #loading the package into the current environment
tcensReg(y ~ 1, data=dt, a=0, v=0.25)
```

    ## $theta
    ##               Estimate
    ## (Intercept)  0.4733706
    ## log_sigma   -0.6936526
    ## 
    ## $iterations
    ## [1] 5
    ## 
    ## $initial_ll
    ## [1] -667.0479
    ## 
    ## $final_ll
    ## [1] -651.3416
    ## 
    ## $var_cov
    ##               (Intercept)     log_sigma
    ## (Intercept)  0.0009230291 -0.0009247406
    ## log_sigma   -0.0009247406  0.0016644569

Note that the this will return parameter estimates, variance-covariance
matrix, the number of iterations until convergence, and the
initial/final log-likelihood values.

Comparing the values to the truth we see that the estimates are
unbiased.

``` r
output <- tcensReg(y ~ 1, data=dt, a=a, v=nu)
lm_output <- lm(y ~ 1, data=dt) #running OLS model for comparison
cens_output <- tcensReg(y ~ 1, data=dt, v=nu) #censored only model, i.e., Tobit model
```

    ## Warning: `a` is not specified indicating no truncation

``` r
tcensReg_est <- output$theta #extracting the point estimates
tcensReg_est[2] <- exp(tcensReg_est[2]) #exponentiating the estimate of log_sigma to estimate sigma

lm_est <- c(coef(lm_output), summary(lm_output)$sigma)

cens_est <- cens_output$theta
cens_est[2] <- exp(cens_est[2])

results_df <- data.frame(rbind(c(mu, sigma), t(tcensReg_est), lm_est, t(cens_est)))
names(results_df) <- c("mu", "sigma")
row.names(results_df) <- c("Truth", "tcensReg", "Normal MLE", "Tobit")
results_df$mu_bias <- abs(results_df$mu - mu)
results_df$sigma_bias <- abs(results_df$sigma - sigma)

knitr::kable(results_df, format="markdown", digits=4)
```

|            |     mu |  sigma | mu\_bias | sigma\_bias |
| :--------- | -----: | -----: | -------: | ----------: |
| Truth      | 0.5000 | 0.5000 |   0.0000 |      0.0000 |
| tcensReg   | 0.4734 | 0.4997 |   0.0266 |      0.0003 |
| Normal MLE | 0.6491 | 0.3645 |   0.1491 |      0.1355 |
| Tobit      | 0.6038 | 0.4307 |   0.1038 |      0.0693 |

Other methods result in significant bias for both `mu` and `sigma`.
