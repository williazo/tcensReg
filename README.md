Maximum Likelihood Estimation of a Truncated Normal Distribution with
Censored Data
================

<!-- badges: start -->

[![CRAN
version](http://www.r-pkg.org/badges/version/tcensReg)](http://www.r-pkg.org/pkg/tcensReg)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/grand-total/tcensReg)](http://www.r-pkg.org/pkg/tcensReg)
[![cran
checks](https://cranchecks.info/badges/worst/tcensReg)](https://cranchecks.info/pkgs/tcensReg)
<!-- badges: end -->

The goal of this package is to estimate parameters from a linear model
when the data comes from a truncated normal distribution with censoring.
Maximum likelihood values are returned. There are multiple method
available for optimization with the default set as conjugate gradient.
This package is also able to return maximum likelihood estimates for
truncated only or censored only data similar to `truncreg` and `censReg`
packages.

For details on the methods, simulation performance, and an application
to visual quality testing from clinical trials for intraocular lenses
see the pre-print manuscript available at
arXiv:[1911.11221](https://arxiv.org/abs/1911.11221).

# Installation

You can install `tcensReg` from CRAN for the stable release or install
the GitHub verion for the active development version:

``` r
#stable CRAN version
install.packages("tcensReg")
#active devel. GitHub version
install.packages("devtools")
devtools::install_github("williazo/tcensReg")
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
mu <- 0.5
sigma <- 0.5
a <- 0
#generate random values from the truncated normal distribution using tcensReg function
y_star <- rtnorm(n=1000, mu=mu, sd=sigma, a=a)
#note that the lowerbound will always be non-negative
round(range(y_star), 3)
```

    ## [1] 0.004 2.355

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

    ## [1] 0.191

``` r
#collecting the uncensored and censored data together
dt <- data.frame(y_star, y) 
```

We can observe the histogram and density plot for the uncensored data,
which shows the zero-truncation.
![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

We can then compare this to the censored observations below
![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

We can then estimate the mean, `mu`, and standard deviation, `sigma`,
using `y` with the `tcensReg` package as shown below.

``` r
#loading the package
library(tcensReg)  
tcensReg(y ~ 1, data=dt, a=0, v=0.25)
```

    ## $theta
    ##               Estimate
    ## (Intercept)  0.5135531
    ## log_sigma   -0.6827764
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $initial_ll
    ## [1] -703.9048
    ## 
    ## $final_ll
    ## [1] -690.939
    ## 
    ## $var_cov
    ##               (Intercept)     log_sigma
    ## (Intercept)  0.0008405050 -0.0008258552
    ## log_sigma   -0.0008258552  0.0015649081
    ## 
    ## $method
    ## [1] "CG"

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
| tcensReg   | 0.5136 | 0.5052 |   0.0136 |      0.0052 |
| Normal MLE | 0.6782 | 0.3740 |   0.1782 |      0.1260 |
| Tobit      | 0.6311 | 0.4445 |   0.1311 |      0.0555 |

Other methods result in significant bias for both `mu` and `sigma`.

# Example 2: Two Population Model with Separate Variances

As an extension for the single population model above, we can imagine a
two independent truncated normal random variables that have common
censoring and truncation values but different standard deviations.

We can simulate the underlying truncated normal distributions `Y1_star`
and `Y2_star` similar to \(Y\) above except now we allow them to have
separate mean and variances.

For this example we let `mu_1`=0.5, `mu_2`=1, `sigma_1`=0.25,
`sigma_2`=2, and `a`=0.

``` r
mu_1 <- 0.5
mu_2 <- 1
sigma_1 <- 0.25
sigma_2 <- 2
a <- 0

y_1_star <- rtnorm(1000, mu = mu_1, sd = sigma_1, a = a)
y_2_star <- rtnorm(1000, mu = mu_2, sd = sigma_2, a = a)
df <- data.frame(y_star = c(y_1_star, y_2_star), 
                 group = c(rep("Population 1", length(y_1_star)),
                           rep("Population 2", length(y_2_star))))
```

Plotting each of these uncensored population densities, we can see the
difference in shape based on the underlying parameter selection.

![](README_files/figure-gfm/two_pop_graph-1.png)<!-- -->

Then censoring each observation at `nu`, we are left with `Y1` and `Y2`.
Again, we let `nu`=0.25.

``` r
nu <- 0.25
df$y <- ifelse(df$y_star<=nu, nu, df$y_star)
```

We then can fit our model with separate variances for each group using
the command `tcensReg_sepvar` as shown below.

``` r
mod_result <- tcensReg_sepvar(y ~ group, a=a, v=nu, group_var="group", method="maxLik", data=df)
mod_result
```

    ## $theta
    ##       (Intercept) groupPopulation 2        log_sigma1        log_sigma2 
    ##        0.50256296       -0.03487324       -1.39156717        0.81610577 
    ## 
    ## $convergence
    ## [1] 2
    ## 
    ## $initial_ll
    ## [1] -2000.264
    ## 
    ## $final_ll
    ## [1] -1900.638
    ## 
    ## $var_cov
    ##                     (Intercept) groupPopulation 2    log_sigma1
    ## (Intercept)        7.628463e-05     -7.628463e-05 -6.926603e-05
    ## groupPopulation 2 -7.628463e-05      7.487459e-02  6.926603e-05
    ## log_sigma1        -6.926603e-05      6.926603e-05  8.017440e-04
    ## log_sigma2        -1.677884e-14     -1.428259e-02  1.523510e-14
    ##                      log_sigma2
    ## (Intercept)       -1.677884e-14
    ## groupPopulation 2 -1.428259e-02
    ## log_sigma1         1.523510e-14
    ## log_sigma2         3.271884e-03
    ## 
    ## $method
    ## [1] "maxLik"

``` r
sepvar_est <- mod_result$theta
sepvar_est[3:4] <- exp(sepvar_est[3:4])

results_df <- data.frame(rbind(c(mu_1, mu_2, sigma_1, sigma_2),
                               t(sepvar_est)))
names(results_df) <- c("mu_1", "mu_2", "sigma_1", "sigma_2")
row.names(results_df) <- c("Truth", "tcensReg")
results_df$mu1_bias <- abs(results_df$mu_1 - mu_1)
results_df$mu2_bias <- abs(results_df$mu_2 - mu_2)
results_df$sigma1_bias <- abs(results_df$sigma_1 - sigma_1)
results_df$sigma2_bias <- abs(results_df$sigma_2 - sigma_2)

knitr::kable(results_df, format="markdown", digits=4)
```

|          |  mu\_1 |    mu\_2 | sigma\_1 | sigma\_2 | mu1\_bias | mu2\_bias | sigma1\_bias | sigma2\_bias |
| :------- | -----: | -------: | -------: | -------: | --------: | --------: | -----------: | -----------: |
| Truth    | 0.5000 |   1.0000 |   0.2500 |   2.0000 |    0.0000 |    0.0000 |       0.0000 |       0.0000 |
| tcensReg | 0.5026 | \-0.0349 |   0.2487 |   2.2617 |    0.0026 |    1.0349 |       0.0013 |       0.2617 |

# Performance Comparison: Censored-Only and Truncated-Only

Note also that the `tcensReg` can also estimate parameters in the
censored-only or truncated-only cases. We show below that by using
analytic values in the tcensReg implementation that our method is faster
then the alternative estimation procedures while providing better
variance estimates. With a small set of covariates and `p<<n` we can use
the Newton Raphson method of optimization, which is computationally fast
with few covariates.

``` r
library(microbenchmark)
#testing the censored-only regression
library(censReg)
cens <- microbenchmark(tcensReg_method = tcensReg(y ~ 1, data=dt, v=nu, method="Newton"),
               censReg_method = censReg(y ~ 1, left=nu, data=dt))
knitr::kable(summary(cens), format="markdown", digits=4)
```

| expr             |     min |      lq |    mean |  median |      uq |     max | neval | cld |
| :--------------- | ------: | ------: | ------: | ------: | ------: | ------: | ----: | :-- |
| tcensReg\_method |  4.8716 |  5.1656 |  7.2273 |  5.6143 |  6.8531 | 22.4004 |   100 | a   |
| censReg\_method  | 13.6368 | 14.6194 | 22.0925 | 17.7750 | 22.5968 | 94.8422 |   100 | b   |

``` r
#point estimates are equivalent
tcensReg_est <- as.numeric(tcensReg(y ~ 1, data=dt, v=nu, method="Newton")$theta)
censReg_est <- as.numeric(coef(censReg(y ~ 1, left=nu, data=dt)))
all.equal(tcensReg_est, censReg_est)
```

    ## [1] TRUE

``` r
#testing the truncated-only regression
library(truncreg)
trunc <- microbenchmark(
  tcensReg_method = tcensReg(y_star ~ 1, data=dt, a=a, method="Newton"),
  truncreg_method = truncreg(y_star ~ 1, point=a, data=dt))
knitr::kable(summary(trunc), format="markdown", digits=4)
```

| expr             |     min |      lq |    mean |  median |      uq |      max | neval | cld |
| :--------------- | ------: | ------: | ------: | ------: | ------: | -------: | ----: | :-- |
| tcensReg\_method |  8.3197 |  9.2489 | 12.8959 | 10.1725 | 14.2182 |  50.7163 |   100 | a   |
| truncreg\_method | 25.2598 | 28.2049 | 36.0450 | 30.8232 | 37.3013 | 168.7519 |   100 | b   |

``` r
tcensReg_est <- as.numeric(tcensReg(y_star ~ 1, data=dt, a=a, method="Newton")$theta)
#note truncreg returns sigma not log_sigma so we need to exponentiate our value
tcensReg_est[2] <- exp(tcensReg_est[2])
truncreg_est <- as.numeric(coef(truncreg(y_star ~ 1, point=a, data=dt)))
all.equal(tcensReg_est, truncreg_est)
```

    ## [1] "Mean relative difference: 1.418798e-07"

In the comparisons above we are using an intercept only model, but in
general we expect that interest lies in understanding how a set of
covariates effect the mean response. So to test the sensitivity and
speed as the number of covariates approaches `n` we can generate
independent random variables `X` and fit the regression model of `Y` or
`Y_star`.

We can compare the censored-only and truncated-only performance with 100
predictors, i.e. `p`=20. To illustrate some of the other available
optimization methods we will set method to BFGS, which is a quasi-Newton
optimization method.

``` r
#number of predictors
p <- 20
X <- NULL
for(i in seq_len(p)){
    X_i <- rnorm(n = length(y))
    X <- cbind(X, X_i)
}
colnames(X) <- paste0("var_", seq_len(p))
dt <- data.frame(y, X)

#testing the censored-only regression with 100 covariates
cens <- microbenchmark(tcensReg_method = tcensReg(y ~ ., data=dt, v=nu, method="BFGS"),
               censReg_method = censReg(y ~ ., left=nu, data=dt))
knitr::kable(summary(cens), format="markdown", digits=4)
```

| expr             |      min |       lq |     mean |   median |       uq |      max | neval | cld |
| :--------------- | -------: | -------: | -------: | -------: | -------: | -------: | ----: | :-- |
| tcensReg\_method | 239.4399 | 259.8322 | 300.1356 | 269.1833 | 323.7974 | 690.5391 |   100 | a   |
| censReg\_method  | 294.4390 | 318.8463 | 375.1485 | 329.4052 | 432.0527 | 665.0162 |   100 | b   |

``` r
#point estimates are equivalent
tcensReg_est <- as.numeric(tcensReg(y ~ ., data=dt, v=nu, method="BFGS")$theta)
censReg_est <- as.numeric(coef(censReg(y ~ ., left=nu, data=dt)))
all.equal(tcensReg_est, censReg_est)
```

    ## [1] "Mean relative difference: 7.673209e-05"

``` r
#testing the truncated-only regression with 100 covariates
trunc <- microbenchmark(tcensReg_method = tcensReg(y_star ~ ., data=dt, a=a, method="BFGS"),
                        truncreg_method = truncreg(y_star ~ ., point=a, data=dt))
knitr::kable(summary(trunc), format="markdown", digits=4)
```

| expr             |      min |       lq |     mean |   median |       uq |      max | neval | cld |
| :--------------- | -------: | -------: | -------: | -------: | -------: | -------: | ----: | :-- |
| tcensReg\_method | 292.8907 | 305.8579 | 333.4667 | 310.8038 | 328.3851 | 693.7038 |   100 | a   |
| truncreg\_method | 339.9185 | 364.1311 | 396.5839 | 370.2025 | 390.6252 | 687.3976 |   100 | b   |

``` r
tcensReg_est <- as.numeric(tcensReg(y_star ~ ., data=dt, a=a, method="BFGS")$theta)
#note truncreg returns sigma not log_sigma so we need to exponentiate the last parameter value
tcensReg_est[length(tcensReg_est)] <- exp(tcensReg_est[length(tcensReg_est)])
truncreg_est <- as.numeric(coef(truncreg(y_star ~ ., point=a, data=dt)))
all.equal(tcensReg_est, truncreg_est)
```

    ## [1] "Mean relative difference: 4.119702e-05"
