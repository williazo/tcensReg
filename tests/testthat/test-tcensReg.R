test_that("error and warning checks", {
    set.seed(1)
    Y <- rnorm(100)
    Y <- ifelse(Y <= 0, 0, Y)
    Y <- ifelse(Y > 1, 1, Y)
    #misordered censoring
    expect_error(tcensReg(Y ~ 1, v = 1, xi = 0),
                   "Left censoring greater than or equal to right censoring")
    #misordered truncation
    expect_error(tcensReg(Y ~ 1, a = 1, b = 0),
                   "Left truncation greater than or equal to right truncation")

    #Censoring specified but no censored values observed
    expect_error(tcensReg(Y ~ 1, v = -1),
                   "Left censoring specified but no censored observations")
    expect_error(tcensReg(Y ~ 1, xi = 2),
                   "Right censoring specified but no censored observations")

    #improper censoring or truncation arguments
    expect_error(tcensReg(Y ~ 1, v = c(0, 2)),
                   "Left censoring value must be numeric scalar")
    expect_error(tcensReg(Y ~ 1, xi = c(0, 2)),
                   "Right censoring value must be numeric scalar")
    expect_error(tcensReg(Y ~ 1, a = "0"),
                   "Left truncation value must be numeric scalar")
    expect_error(tcensReg(Y ~ 1, b = list(0, 1)),
                   "Right truncation value must be numeric scalar")

    #no censoring or truncation warning
    expect_message(tcensReg(Y ~ 1),
                   "No censoring and no truncation")

    #left/right truncation
    #improper truncation value
    expect_error(tcensReg(Y ~ 1, a = -1, v = 0, xi = 1, b = 0.5),
                   "Observed values outside truncation value")
    expect_error(tcensReg(Y ~ 1, a = 0.5, v = 0, xi = 1, b = 2),
                   "Observed values outside truncation value")
    #proper specification
    expect_message(tcensReg(formula = Y ~ 1, a = -1, v = 0, xi = 1, b = 2),
                   "Left/right censoring with left/right truncation")
    #improper order
    expect_error(tcensReg(Y ~ 1, a = 0, v = 0, xi = 1, b = 2),
                   "Censoring specified after truncation value")
    #improper order
    expect_error(tcensReg(Y ~ 1, a = 0, v = 0, xi = 1, b = 1),
                   "Censoring specified after truncation value")

    #censoring only, no truncation
    expect_message(tcensReg(Y ~ 1, v = 0),
                    "Left censoring with no truncation")
    expect_message(tcensReg(Y ~ 1, v = 0, xi = 1),
                    "Left/right censoring with no truncation")
    expect_message(tcensReg(Y ~ 1, xi = 1),
                    "Right censoring with no truncation")

    #left truncation
    expect_error(tcensReg(Y ~ 1, a = 0.5),
                    "Observed values below left truncation value")
    expect_error(tcensReg(Y ~ 1, a = 0, v = 0, xi = 1),
                    "Censoring specified below/equal to left truncation")
    expect_warning(tcensReg(Y ~ 1, a = 0, xi = 1),
                    "Right censoring with left truncation")

    #right truncation
    expect_error(tcensReg(Y ~ 1, b = 0.5),
                    "Observed values above right truncation value")
    expect_error(tcensReg(Y ~ 1, v = 0, xi = 1, b = 1),
                    "Censoring specified above/equal to right truncation")
    expect_warning(tcensReg(Y ~ 1, v = 0, b = 2),
                    "Left censoring with right truncation")
})

test_that("functional tcensReg output check", {
    y_star <- rtnorm(n = 1000, mu = 0.5, sd = 1, a = 0)
    y <- ifelse(y_star <= 0.25, 0.25, y_star)
    output <- tcensReg(y ~ 1, v=0.25, a=0)

    #correct number of objects
    expect_equal(length(output), 11)
    #correct dimensions of variance covariance
    expect_equal(dim(output$var_cov), c(2, 2))
    #final log-likelihood above initial
    expect_true(output$final_ll >  output$initial_ll)
})
