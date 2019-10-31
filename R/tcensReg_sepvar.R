#' @title Regression Method for Truncated Normal Distribution with Censoring for Indepentend Truncated Normal Groups
#'
#' @description This function is used to find estimates from a linear equation assuming that the data is observed from a truncated
#'  distribution with left censoring. It uses numerical values of the gradient vector and hessian matrix to
#'  solve for the maximum likelihood using \code{maxLik} package. This function can also
#'  be used with censored only, truncated only, or uncensored and untruncated gaussian models.
#'
#' @param formula Object of class \code{formula} which symbolically describes the model to be fit
#' @param a Numeric scalar indicating the truncation value. Initial value is -Inf indicating no truncation
#' @param v Numeric scalar indicating the censoring value. Initially set to NULL indicating no censoring
#' @param group_var Character scalar indiciating a variable in the data.frame that defines the independent groups
#' @param method Character value indicating which optimization routine to perform. Choices include \code{BFGS}, \code{maxLik} and \code{CG}. See details for explanation on each method.
#' @param theta_init Optional initial values for the parameters. Default is to fit a linear regression model.
#' @param data Data.frame that contains the outcome and corresponding covariates. If none is provided then assumes objects are in user's environment.
#' @param max_iter Numeric value indicating the maximum number of iterations to perform.
#' @param ... Additional arguments from \code{\link{tcensReg_newton}} such as \code{max_iter}, \code{step_max}, or \code{epsilon}.
#'
#' @importFrom stats model.frame model.matrix
#' @importFrom maxLik maxLik
#'
#' @details
#'  Currently available optimization routines include conjugate gradient (\code{CG}), Newton-Raphson type via maxLike package (\code{\link[maxLik]{maxLik}}), and BFGS (\code{BFGS}).
#'  The default method is set as the conjugate gradient. Both the of the conjugate gradient and BFGS methods are implemented via the
#'  general-purpose optimization \code{\link{optim}}. These two methods use only the respective likelihood and gradient functions.
#'  The Newton-Raphson method uses the likelihood, gradient, and Hessian functions along with line search to achieve the maximum likelihood.
#'
#' @return Returns a list of final estimate of theta, total number of iterations performed, initial log-likelihood,
#' final log-likelihood, and estimated variance covariance matrix.
#'
#' @export

tcensReg_sepvar <- function(formula, a = -Inf, v = NULL, group_var,
                            method=c("BFGS", "maxLik", "CG"),
                            theta_init = NULL, data = sys.frame(sys.parent()),
                            max_iter = 100,
                            ...){
  method <- match.arg(method)
  #checks for proper specification of formula
  if(class(formula)!= "formula"){
    stop("`formula` must be a formula", call. = FALSE)
  } else if (length(formula) != 3){
    stop("`formula` must be 2-sided", call. = FALSE)
  }

  y <- model.frame(formula, data)[, 1]
  group <- data[, group_var]
  if(!is.factor(group)) group <- as.factor(as.character(group))
  X <- model.matrix(formula, data)


  #checking for proper specification of a and v
  if(!is.null(v) & (!is.numeric(a) | !is.numeric(v))){
    stop("`a` and `v` must both be numeric", call. = FALSE)
  } else if(!is.null(v) & (length(a) !=1 | length(v) != 1)){
    stop("`a`, and `v` must both be scalars", call. = FALSE)
  }
  # checking model specification of truncation and censoring params
  if(is.null(v) & a == -Inf){
    warning("`v` and `a` are not specified indicating no censoring and no truncation", call. = FALSE)
  } else if(is.null(v) & a != -Inf){
    if(any(y < a)){
      stop("observed values below specified truncation `a`", call.=FALSE)
    }
    warning("`v`is not specified indicating no censoring", call. = FALSE)
  } else if(!is.null(v) & a == -Inf){
    len_cens <- sum(y == v, na.rm=TRUE)
    if(len_cens == 0){
      stop("censoring indicated but no observed censored values",
           call.=FALSE)
    }
    warning("`a` is not specified indicating no truncation", call. = FALSE)
  } else {
    if(v < a) stop("censoring specified below truncation", call.=FALSE)
    if(any(y < a)){
      stop("observed values below specified truncation `a`", call.=FALSE)
    }
    len_cens <- sum(y == v, na.rm = TRUE)
    if(len_cens == 0){
      stop("censoring indicated but no observed censored values",
           call.=FALSE)
    }
  }


  if(!group_var%in%names(model.frame(data))){
    stop("`group_var` must be a variable in the data.frame", call. = FALSE)
  }

  #want to use different inital estimates depending on whether it is truncation only, censor only, or truncated and censored
  #if censored only, normal, or truncated only then use estimates from OLS
  if(is.null(theta_init)==TRUE){
    lm_mod <- lm(y ~ X - 1)
    beta_init <- unname(coef(lm_mod))

    num_groups <- length(unique(group))
    log_sigmas <- vector(length = num_groups)
    for(j in 1:num_groups){
      log_sigmas[j] <- log(sd(y[group == unique(group)[j]]))
    }
    theta_init <- c(beta_init, log_sigmas)
  }

  names(theta_init)[1:(length(theta_init)-length(unique(group)))] <- colnames(X)
  names(theta_init)[(length(theta_init)-length(unique(group))+1):length(theta_init)] <- paste0("log_sigma", 1:length(unique(group)))

  #reading in the newton raphson for the truncated censored normal
  if(method %in% c("BFGS", "CG")){
    null_ll <- tcensReg_llike_sepvar(theta_init, y, X, group, a, v)

    optim_results <- optim(par=theta_init,
                           fn=tcensReg_llike_sepvar, gr=tcensReg_gradient_sepvar,
                           method=method, hessian=TRUE,
                           y=y, X=X, a=a, v=v, group=group,
                           control=list(maxit = max_iter, fnscale=-1))
    if(optim_results$convergence==1){
        warning("Convergence unsuccessful: maximum iterations reached", call.=FALSE)
      }
    v_cov <- -solve(optim_results$hessian) #converting the hessian to estimate of variance
    theta <- matrix(optim_results$par, nrow=length(optim_results$par), ncol=1)

    # row.names(theta) <- c(colnames(X),"log_sigma") #adding in variable names
    colnames(theta) <- "Estimate"
    # row.names(v_cov) <- row.names(theta) #using the names for the variance covariance matrix
    # colnames(v_cov) <- row.names(theta)
    results <- list(theta = theta, convergence = optim_results$convergence,
                initial_ll = null_ll, final_ll = optim_results$value, var_cov = v_cov,
                method=method)
  }
  else if(method == "maxLik"){
    left_trunc <- a
    null_ll <- tcensReg_llike_sepvar_maxLik(theta_init, y, X, group, a, v)

    maxLik_results <- maxLik::maxLik(logLik=tcensReg_llike_sepvar_maxLik,
                                     grad=tcensReg_gradient_sepvar_maxLik,
                                     start=theta_init,
                                     y=y, X=X, left_trunc=left_trunc, v=v, group=group)
    results <- list(theta = coef(maxLik_results), convergence = maxLik_results$code,
                    initial_ll = null_ll, final_ll=summary(maxLik_results)$loglik,
                    var_cov = vcov(maxLik_results), method = "maxLik")

  }

  return(results)
}
