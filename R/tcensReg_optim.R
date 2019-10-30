#' Optimization for Truncated Normal Distribution with Censoring with Linear Equation Mean using optim
#'
#' @param y Numeric vector with the observed truncated and censored outcomes
#' @param X Numeric design matrix
#' @param a Numeric scalar indicating the truncation value. Initial value is -Inf indicating no truncation
#' @param v Numeric scalar indicating the censoring value. Initially set to NULL indicating no censoring
#' @param method Character value indicating which optimization method to perform.
#' @param epsilon Numeric value used to define when the algorithm should stop when the gradient is less then epsilon. Default is 0.001
#' @param theta_init Initial values of theta provided by the user. If unspecified then calculcates values from OLS regression
#' @param max_iter Maximum number of iterations for algorithm. Default is 100
#' @param step_max Maximum number of steps when performing line search. Default is 10
#' @param tol_val Tolerance value used to stop the algorithm if the (n+1) and (n) log likelihood is within the tolerance limit
#'
#' @importFrom stats coef dnorm lm model.frame model.matrix pnorm optim
#'
#' @export
#'
#' @return Returns a list of final estimate of theta, total number of iterations performed, initial log-likelihood,
#' final log-likelihood, and estimated variance covariance matrix.


tcensReg_optim<-function(y, X, a = -Inf, v = NULL, method, epsilon = 1e-4,
                          tol_val = 1e-6, max_iter = 100, step_max = 10,
                          theta_init = NULL){

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
