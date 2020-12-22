#' Gradient Vector for Truncated Normal Distribution with Censoring with Linear Equation Mean
#'
#' @inheritParams tcensReg_llike
#'
#' @importFrom stats dnorm pnorm
#' @return Vector of gradient values with p-1 beta parameters and log sigma for the nth iterate
#' @keywords internal
tcensReg_gradient <- function(theta,
                              y,
                              X,
                              a = -Inf,
                              v = NULL,
                              xi = NULL,
                              b = Inf) {

  nabla <- vector(length = length(theta)) #creating an empty vector to store the gradient

  #assume that there are a total of p parameters, then the first p-1 are assumed to be beta
  #and the last parameter is log_sigma
  p <- length(theta)

  log_sigma <- theta[p]

  n <- length(y) #total number of observations

  if(is.null(v) == FALSE){
    uncens <- which(y > v) #identifying which values are uncensored
    cens <- which(y == v) #idenfitying which values are censored

    n_0 <- length(cens) #number of censored observations
    n_1 <- length(uncens) #number of uncensored observations
    a_stand <- (a - X %*% theta[seq_len(p - 1)]) / exp(log_sigma) #standardized value with respect to truncated value
    v_stand <- (v - X %*% theta[seq_len(p - 1)]) / exp(log_sigma) #standardized value with respect to censored value
    y_stand <- (y - X %*% theta[seq_len(p - 1)]) / exp(log_sigma)

    #calculating the gradient for each beta parameter
    nabla_beta <- vapply(seq_len(p - 1), FUN.VALUE = numeric(1), FUN = function(j) {
      d1 <- sum((-X[, j] * exp(-log_sigma) * dnorm(a_stand)) / pnorm(-a_stand))
      d2 <- sum((X[cens, j] * exp(-log_sigma) * (dnorm(v_stand[cens]) - dnorm(a_stand[cens]))) / (pnorm(v_stand[cens]) - pnorm(a_stand[cens])))
      d3 <- sum(X[uncens, j] * exp(-log_sigma) * (y_stand[uncens]))
      nabla <- d1 - d2 + d3
      return(nabla)
    })

    #here we need a special case when there is no truncation (i.e. a = -Inf)
    if (a == -Inf){
      d1 <- 0
      d2 <- sum((v_stand[cens] * dnorm(v_stand[cens])) / (pnorm(v_stand[cens])))
    } else{
      d1 <- sum((-a_stand * dnorm(a_stand)) / pnorm(-a_stand))
      d2 <- sum((v_stand[cens] * dnorm(v_stand[cens])-a_stand[cens] * dnorm(a_stand[cens])) / (pnorm(v_stand[cens]) - pnorm(a_stand[cens])))
    }
    d3 <- n_1
    d4 <- sum(y_stand[uncens] ^ 2)
    nabla_sigma <- d1 - d2 - d3 + d4

    nabla <- c(nabla_beta, nabla_sigma)
  } else if (is.null(v) == TRUE){
    a_stand <- (a - X %*% theta[seq_len(p - 1)]) / exp(log_sigma) #standardized value with respect to truncated value
    y_stand <- (y - X %*% theta[seq_len(p - 1)]) / exp(log_sigma)

    #calculating the gradient for each beta parameter
    nabla_beta <- vapply(seq_len(p - 1), FUN.VALUE = numeric(1), FUN = function(j){
      d1 <- sum((-X[, j] * exp(-log_sigma) * dnorm(a_stand)) / pnorm(-a_stand))
      d3 <- sum(X[, j] * exp(-log_sigma) * (y_stand))
      nabla <- d1 + d3
      return(nabla)
    })

    if(a == -Inf){
      d1 <- 0
    }else{
      d1 <- sum((-a_stand * dnorm(a_stand)) / pnorm(-a_stand))
    }
    d3 <- n
    d4 <- sum(y_stand ^ 2)
    nabla_sigma <- d1 - d3 + d4
    nabla <- c(nabla_beta, nabla_sigma)
  }
  return(nabla)
}
