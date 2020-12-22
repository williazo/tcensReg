#' Optimization for Truncated Normal Distribution with Censoring with Linear Equation Mean using optim
#'
#' @description Iteratively solves likelihood of truncated normal with censoring
#' using only gradient and/or likelihood values
#'
#' @inheritParams tcensReg_llike
#' @inheritParams tcensReg
#' @inheritParams tcensReg_newton
#'
#' @importFrom stats coef dnorm lm model.frame model.matrix pnorm optim
#'
#' @return Returns a list of final estimate of theta, total number of iterations performed, initial log-likelihood,
#' final log-likelihood, and estimated variance covariance matrix.
#' @keywords internal


tcensReg_optim<-function(y,
                         X,
                         a = -Inf,
                         v = NULL,
                         xi = NULL,
                         b = Inf,
                         method,
                         epsilon = 1e-4,
                         tol_val = 1e-6,
                         max_iter = 100,
                         step_max = 10,
                         theta_init = NULL) {

    #want to use different inital estimates depending on whether it is truncation only, censor only, or truncated and censored
    #if censored only, normal, or truncated only then use estimates from OLS
    if(is.null(theta_init) == TRUE & (a != -Inf & is.null(v) == FALSE)){
        #if censored and truncated then use initial estimates from censored only model
        cens_mod <- suppressWarnings(tcensReg(y ~ X - 1, v = v, method=method))
        theta_init <- unname(cens_mod$theta)
    } else{
        lm_mod <- lm(y ~ X - 1)
        theta_init <- c(unname(coef(lm_mod)), log(unname(summary(lm_mod)$sigma)))
    }

    null_ll <- tcensReg_llike(theta_init, y, X, a, v)

    optim_results <- optim(par=theta_init,
                          fn=tcensReg_llike, gr=tcensReg_gradient,
                          method=method, hessian=TRUE,
                          y=y, X=X, a=a, v=v,
                          control=list(maxit=max_iter, fnscale=-1))
    v_cov <- -solve(optim_results$hessian) #converting the hessian to estimate of variance
    theta <- matrix(optim_results$par, nrow=length(optim_results$par), ncol=1)

    row.names(theta) <- c(colnames(X),"log_sigma") #adding in variable names
    colnames(theta) <- "Estimate"
    row.names(v_cov) <- row.names(theta) #using the names for the variance covariance matrix
    colnames(v_cov) <- row.names(theta)

    return(list(theta = theta, convergence = optim_results$convergence,
                initial_ll = null_ll, final_ll = optim_results$value, var_cov = v_cov,
                method=method))
}
