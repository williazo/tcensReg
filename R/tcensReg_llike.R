#' Log Likelihood for Truncated Normal Distribution with Censoring with Linear Equation Mean
#'
#' \code{tcensReg_llike} is a supporting function that is used to calculate the log likelihood value at the nth
#' iteration of theta. If \code{a} and/or \code{v} are not specified then the corresponding censored only, truncated only,
#' or gaussian log likelihood will be used. This function is called as part of the Newton-Raphson algorithm in \code{tcensReg_newton}.
#'
#' @param theta Numeric vector numeric vector containing estimates of beta and log sigma
#' @param a Numeric scalar indicating the truncation value
#' @param v Numeric scalar indicating the censoring value
#' @param y Numeric vector with the observed truncated and censored outcomes
#' @param X Numeric design matrix
#'
#' @importFrom stats dnorm pnorm
#' @export
#'
#' @return Scalar value of the log-likelihood at the nth iterate

tcensReg_llike <- function(theta, y, X, a = -Inf, v = NULL){

  #assume that there are a total of p parameters, then the first p-1 are assumed to be beta
  #and the last parameter is log_sigma
  # if(min(y)<min(a, v)){stop("Values found below the censored or truncation value", call. = FALSE)}
  p <- length(theta)
  log_sigma <- theta[p]

  n <- length(y) #total number of observations

  if(is.null(v)==FALSE){
    uncens <- which(y>v) #identifying which values are uncensored
    cens <- which(y == v) #idenfitying which values are censored
    if(min(y)<min(a, v)){stop("Values found below the censored or truncation value", call. = FALSE)}
    if(v<=a){stop("Censoring value cannot be below truncation value", call. = FALSE)}
    n_0 <- length(cens) #number of censored observations
    n_1 <- length(uncens) #number of uncensored observations
    a_stand <- (a-X%*%theta[1:(p-1)])/exp(log_sigma) #standardized value with respect to truncated value
    v_stand <- (v-X%*%theta[1:(p-1)])/exp(log_sigma) #standardized value with respect to censored value
    y_stand <- (y-X%*%theta[1:(p-1)])/exp(log_sigma)

    l_lik <- -sum(log(pnorm(-a_stand))) + sum(log(pnorm(v_stand[cens])-pnorm(a_stand[cens]))) - n_1*log_sigma + sum(log(dnorm(y_stand[uncens])))
  }else if(is.null(v)==TRUE){
    a_stand <- (a-X%*%theta[1:(p-1)])/exp(log_sigma) #standardized value with respect to truncated value
    y_stand <- (y-X%*%theta[1:(p-1)])/exp(log_sigma)

    l_lik <- -sum(log(pnorm(-a_stand))) - n*log_sigma + sum(log(dnorm(y_stand)))
  }
  return(l_lik)
}
