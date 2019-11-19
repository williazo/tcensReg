test_that("error and warning checks", {
    Y <- rnorm(100)
    Y <- ifelse(Y <= 0, 0, Y)
    #truncation only warning
    expect_warning(tcensReg(Y ~ 1, v = 0),
                    "`a` is not specified indicating no truncation")
    #censored only warning
    expect_warning(tcensReg(Y ~ 1, a = -1),
                    "`v`is not specified indicating no censoring")
    #no censoring or truncation warning
    expect_warning(tcensReg(Y ~ 1),
                   "`v` and `a` are not specified indicating no censoring and no truncation")
    #no observed censored values in censored only
    expect_error(tcensReg(Y ~ 1, v=1),
                   "censoring indicated but no observed censored values")
    #no obseved censored values in tcensReg
    expect_error(tcensReg(Y ~ 1, v=1, a=-1),
                 "censoring indicated but no observed censored values")
    #values below truncation observed, truncation only model
    expect_error(tcensReg(Y ~ 1, a=1),
                 "observed values below specified truncation `a`")
    #values below truncation observed, tcensReg model
    expect_error(tcensReg(Y ~ 1, v=10, a=1),
                 "observed values below specified truncation `a`")
    #censoring below trunction
    expect_error(tcensReg(Y ~ 1, v=0, a=1),
                 "censoring specified below truncation")
})

test_that("functional tcensReg output check", {
    y_star <- rtnorm(n = 1000, mu = 0.5, sd = 1, a = 0)
    y <- ifelse(y_star <= 0.25, 0.25, y_star)
    output <- tcensReg(y ~ 1, v=0.25, a=0)

    #correct number of objects
    expect_equal(length(output), 6)
    #correct dimentsions of variance covariance
    expect_equal(dim(output$var_cov), c(2, 2))
    #final log-likelihood below initial
    expect_true(output$final_ll >  output$initial_ll)
})
