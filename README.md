Maximum Likelihood Estimation of a Truncated Normal Distribution with Censored Data
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
devtools::install_github("williazo/tcensReg",
                         build_opts=c("--no-resave-data", "--no-manual"),
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
it as `y_star`.

``` r
library(msm)
mu <- 0.5
sigma <- 0.5
a <- 0
#generate random values from the truncated normal distribution
y_star <- msm::rtnorm(n=1000, mean=mu, sd=sigma, lower=a)
#note that the lowerbound will always be non-negative
round(range(y_star), 3)
```

    ## [1] 0.001 2.026

Next, we can imagine a scenario where we have an imprecise measurement
of `y_star` leading to censoring. In our case we assume that values
below a limit of detection, `nu`, are censored. This creates a random
variable `y`.

In the example below we set our limit of detection as `nu`=0.25.

``` r
nu <- 0.25
y <- ifelse(y_star <= nu, nu, y_star)
#calculating the number of censored observations
sum(y == nu)/length(y)
```

    ## [1] 0.171

``` r
#collecting the uncensored and censored data together
dt <- data.frame(y_star, y)
```

We can observe the histogram and density plot for the uncensored data,
which shows the zero-truncation.
![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

We can then compare this to the censored observations below
![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

We can then estimate the mean, `mu`, and standard deviation, `sigma`,
using `y` with the `tcensReg` package as shown below.

``` r
#loading the package
library(tcensReg)  
tcensReg(y ~ 1, data=dt, a=0, v=0.25)
```

    ## $theta
    ##               Estimate
    ## (Intercept)  0.5435955
    ## log_sigma   -0.7216200
    ##
    ## $iterations
    ## [1] 5
    ##
    ## $initial_ll
    ## [1] -668.4731
    ##
    ## $final_ll
    ## [1] -657.1288
    ##
    ## $var_cov
    ##               (Intercept)     log_sigma
    ## (Intercept)  0.0006611558 -0.0006590367
    ## log_sigma   -0.0006590367  0.0014293139

Note that the this will return parameter estimates, variance-covariance
matrix, the number of iterations until convergence, and the
initial/final log-likelihood values.

Comparing the values to the truth we see that the estimates are
unbiased.

``` r
#tcensReg model
output <- tcensReg(y ~ 1, data=dt, a=a, v=nu)
#extracting the point estimates
tcensReg_est <- output$theta
#exponentiating the estimate of log_sigma to obtain sigma
tcensReg_est[2] <- exp(tcensReg_est[2])

#OLS model
lm_output <- lm(y ~ 1, data=dt)
lm_est <- c(coef(lm_output), summary(lm_output)$sigma)
#censored only model, i.e., Tobit model
cens_output <- tcensReg(y ~ 1, data=dt, v=nu)
```

    ## Warning: `a` is not specified indicating no truncation

``` r
cens_est <- cens_output$theta
cens_est[2] <- exp(cens_est[2])


results_df <- data.frame(rbind(c(mu, sigma),
                               t(tcensReg_est),
                               lm_est,
                               t(cens_est)))
names(results_df) <- c("mu", "sigma")
row.names(results_df) <- c("Truth", "tcensReg", "Normal MLE", "Tobit")
results_df$mu_bias <- abs(results_df$mu - mu)
results_df$sigma_bias <- abs(results_df$sigma - sigma)

knitr::kable(results_df, format="markdown", digits=4)
```

|            |     mu |  sigma | mu\_bias | sigma\_bias |
| :--------- | -----: | -----: | -------: | ----------: |
| Truth      | 0.5000 | 0.5000 |   0.0000 |      0.0000 |
| tcensReg   | 0.5436 | 0.4860 |   0.0436 |      0.0140 |
| Normal MLE | 0.6828 | 0.3704 |   0.1828 |      0.1296 |
| Tobit      | 0.6427 | 0.4316 |   0.1427 |      0.0684 |

Other methods result in significant bias for both `mu` and `sigma`.

Note also that the `tcensReg` can also estimate parameters in the
censored-only or truncated-only cases. We show below that by using
analytic values in the tcensReg implementation that our method is faster
then the alternative estimation procedures while providing better
variance estimates.

``` r
library(microbenchmark)
#testing the censored-only regression
library(censReg)
cens <- microbenchmark(tcensReg_method = tcensReg(y ~ 1, data=dt, v=nu),
               censReg_method = censReg(y ~ 1, left=nu, data=dt))
knitr::kable(summary(cens), format="markdown", digits=4)
```

| expr             |     min |      lq |    mean |  median |      uq |      max | neval |
| :--------------- | ------: | ------: | ------: | ------: | ------: | -------: | ----: |
| tcensReg\_method |  4.7981 |  5.6043 |  7.5102 |  6.4579 |  8.0752 |  28.8798 |   100 |
| censReg\_method  | 15.0433 | 17.9193 | 27.5966 | 22.8280 | 28.0810 | 187.2731 |   100 |

``` r
#point estimates are equivalent
tcensReg_est <- as.numeric(tcensReg(y ~ 1, data=dt, v=nu)$theta)
censReg_est <- as.numeric(coef(censReg(y ~ 1, left=nu, data=dt)))
all.equal(tcensReg_est, censReg_est)
```

    ## [1] TRUE

``` r
#testing the truncated-only regression
library(truncreg)
trunc <- microbenchmark(tcensReg_method = tcensReg(y_star ~ 1, data=dt, a=a),
                        truncreg_method = truncreg(y_star ~ 1, point=a, data=dt))
knitr::kable(summary(trunc), format="markdown", digits=4)
```

| expr             |     min |      lq |    mean |  median |      uq |      max | neval |
| :--------------- | ------: | ------: | ------: | ------: | ------: | -------: | ----: |
| tcensReg\_method |  8.5226 |  9.7616 | 11.9713 | 10.8375 | 13.4695 |  36.5793 |   100 |
| truncreg\_method | 27.3353 | 30.5580 | 38.8867 | 34.0826 | 40.0806 | 175.2784 |   100 |

``` r
tcensReg_est <- as.numeric(tcensReg(y_star ~ 1, data=dt, a=a)$theta)
#note truncreg returns sigma not log_sigma so we need to exponentiate our value
tcensReg_est[2] <- exp(tcensReg_est[2])
truncreg_est <- as.numeric(coef(truncreg(y_star ~ 1, point=a, data=dt)))
all.equal(tcensReg_est, truncreg_est)
```

    ## [1] "Mean relative difference: 2.851958e-08"
