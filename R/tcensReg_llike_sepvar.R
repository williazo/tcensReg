#' Log Likelihood for Model for J independent Truncated Normal Variables with Same Mean Structure but Different Variance
#'
#' This function is part of the supporting function space to calculate the value of the log-likelihood at the nth iteration of theta.
#' For this parametrization we assume that there exists (p-1) linear parameters for the mean and J separate variance parameters.
#'
#' @param theta Numeric vector containing the values at the present iterate for the (p-1) fixed mean parameters and the J log sigma values
#' @param a Numeric scalar defining the common truncation value for each Truncated Normal
#' @param v Numeric scalar defining the common censoring value
#' @param y Numeric vector containing the outcome
#' @param X Numeric design matrix
#' @param group Factor variable used to define the J groups
#'
#' @importFrom stats dnorm pnorm
#' @export
#'
#' @return Scalar value of the log-likelihood at the nth-iterate
tcensReg_llike_sepvar <- function(theta, y, X, group, a = -Inf, v = NULL){
  p <- length(theta)
  num_groups <- length(unique(group))
  log_sigmas <- theta[(p-num_groups+1):p]
  lin_pred <- theta[1:(p-num_groups)]

  loglik_components <- rep(NA, num_groups)
  for(j in 1:num_groups){
    y_j <- y[group == unique(group)[j]]
    X_j <- as.matrix(X[group == unique(group)[j], ], nrow = length(group == unique(group)[j]), ncol = ncol(X)) #I want to ensure this is still a matrix in order to use the matrix multiplicaiton later
    n_j <- length(y_j) #number of observations in jth group

    if(is.null(v)==FALSE){
      uncens <- which(y_j>v) #identifying which values are uncensored
      cens <- which(y_j == v) #idenfitying which values are censored
      if(min(y_j)<min(a, v)){stop("Values found below the censored or truncation value", call. = FALSE)}
      if(v<=a){stop("Censoring value cannot be below truncation value", call. = FALSE)}
      n_0j <- length(cens) #number of censored observations
      n_1j <- length(uncens) #number of uncensored observations

      a_stand <- (a-X_j%*%lin_pred)/exp(log_sigmas[j]) #standardized value with respect to truncated value
      v_stand <- (v-X_j%*%lin_pred)/exp(log_sigmas[j]) #standardized value with respect to censored value
      y_stand <- (y_j-X_j%*%lin_pred)/exp(log_sigmas[j])

      l_lik <- -sum(log(pnorm(-a_stand))) + sum(log(pnorm(v_stand[cens])-pnorm(a_stand[cens]))) - n_1j*log_sigmas[j] + sum(log(dnorm(y_stand[uncens])))
      loglik_components[j] <- l_lik
      }else if(is.null(v)==TRUE){
      a_stand <- (a-X_j%*%lin_pred)/exp(log_sigmas[j]) #standardized value with respect to truncated value
      y_stand <- (y_j-X_j%*%lin_pred)/exp(log_sigma[j])

      l_lik <- -sum(log(pnorm(-a_stand))) - n_j*log_sigma + sum(log(dnorm(y_stand)))
      loglik_components[j] <- l_lik
    }
  }
  loglik <- sum(loglik_components)
  return(loglik)
}
