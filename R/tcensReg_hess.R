#' Hessian Matrix for Truncated Normal Distribution with Censoring with Linear Equation Mean
#' @param theta Numeric vector numeric vector containing estimates of beta and log sigma
#' @param a Numeric scalar indicating the truncation value
#' @param v Numeric scalar indicating the censoring value
#' @param y Numeric vector with the observed truncated and censored outcomes
#' @param X Numeric design matrix
#'
#' @importFrom stats dnorm pnorm
#' @export
#'
#' @return Matrix of Hessian values for the nth iterate

tcensReg_hess <- function(theta, y, X, a = -Inf, v = NULL){

  H <- matrix(nrow = length(theta), ncol = length(theta)) #creating an empty vector to store the gradient

  #assume that there are a total of p parameters, then the first p-1 are assumed to be beta
  #and the last parameter is log_sigma
  p <- length(theta)

  log_sigma <- theta[p]

  n <- length(y) #total number of observations

  if(is.null(v)==FALSE){
    uncens <- which(y>v) #identifying which values are uncensored
    cens <- which(y == v) #idenfitying which values are censored

    n_0 <- length(cens) #number of censored observations
    n_1 <- length(uncens) #number of uncensored observations
    a_stand <- (a-X%*%theta[1:(p-1)])/exp(log_sigma) #standardized value with respect to truncated value
    v_stand <- (v-X%*%theta[1:(p-1)])/exp(log_sigma) #standardized value with respect to censored value
    y_stand <- (y-X%*%theta[1:(p-1)])/exp(log_sigma)


    if(a == -Inf){
      #calculating the hessian for each beta parameter
      for(j in 1:(p-1)){
        for(k in j:(p-1)){
          H1 <- 0
          H2 <- sum((X[cens,j]*X[cens,k]*exp(-2*log_sigma))*(((pnorm(v_stand[cens]))*(dnorm(v_stand[cens])*v_stand[cens])+(dnorm(v_stand[cens]))^2)/(pnorm(v_stand[cens]))^2))
          H3 <- (X[uncens,j]%*%X[uncens,k])*exp(-2*log_sigma)
          H[j, k] <- H1 - H2 - H3
        }
      }
      #calculating the hessian for beta with respect to log_sigma
      for(j in 1:(p-1)){
        H1 <- 0
        H2 <- sum((-X[cens,j]*exp(-log_sigma))*((pnorm(v_stand[cens]))*((dnorm(v_stand[cens])*(1-v_stand[cens]^2)))-((dnorm(v_stand[cens]))*((dnorm(v_stand[cens])*v_stand[cens]))))/((pnorm(v_stand[cens]))^2))
        H3 <- -2*exp(-log_sigma)*(X[uncens,j]%*%y_stand[uncens])
        H[j, p] <- H1 - H2 + H3
      }
      #this completes the upper triangle of the hessian matrix

      #calculating the hessian with respect to log_sigma

      #calculating the second derivative with respect to log_sigma squared
      H1 <- 0
      H2 <- sum(((pnorm(v_stand[cens]))*(((dnorm(v_stand[cens])*(v_stand[cens]^3-v_stand[cens]))))+((dnorm(v_stand[cens])*v_stand[cens]))^2)/((pnorm(v_stand[cens]))^2))
      H3 <- 2*sum((y_stand[uncens]^2))
      H[p,p] <- H1 - H2 - H3

      #completing the matrix using the symmetric property
      H[lower.tri(H)]<- H[upper.tri(H)]
    }else{
      #calculating the hessian with respect to the beta parameters
      for(j in 1:(p-1)){
        for(k in j:(p-1)){
          H1 <- sum((-X[,j]*X[,k]*exp(-2*log_sigma))*((pnorm(-a_stand)*dnorm(a_stand)*a_stand-dnorm(a_stand)^2)/pnorm(-a_stand)^2))
          H2 <- sum((X[cens,j]*X[cens,k]*exp(-2*log_sigma))*(((pnorm(v_stand[cens])-pnorm(a_stand[cens]))*(dnorm(v_stand[cens])*v_stand[cens]-dnorm(a_stand[cens])*a_stand[cens])+(dnorm(v_stand[cens])-dnorm(a_stand[cens]))^2)/(pnorm(v_stand[cens])-pnorm(a_stand[cens]))^2))
          H3 <- (X[uncens,j]%*%X[uncens,k])*exp(-2*log_sigma)
          H[j, k] <- H1 - H2 - H3
        }
      }

      #calculating the hessian for beta with respect to log_sigma
      for(j in 1:(p-1)){
        H1 <- sum((X[,j]*exp(-log_sigma))*((pnorm(-a_stand)*(dnorm(a_stand)*(1-a_stand^2))+a_stand*dnorm(a_stand)^2)/pnorm(-a_stand)^2))
        H2 <- sum((-X[cens,j]*exp(-log_sigma))*((pnorm(v_stand[cens])-pnorm(a_stand[cens]))*((dnorm(v_stand[cens])*(1-v_stand[cens]^2))-(dnorm(a_stand[cens])*(1-a_stand[cens]^2)))-((dnorm(v_stand[cens])-dnorm(a_stand[cens]))*((dnorm(v_stand[cens])*v_stand[cens])-(dnorm(a_stand[cens])*a_stand[cens]))))/((pnorm(v_stand[cens])-pnorm(a_stand[cens]))^2))
        H3 <- -2*exp(-log_sigma)*(X[uncens,j]%*%y_stand[uncens])
        H[j, p] <- H1 - H2 + H3
      }
      #this completes the upper triangle of the hessian matrix

      #calculating the hessian with respect to log_sigma

      #calculating the second derivative with respect to log_sigma squared
      H1 <- sum((a_stand*(pnorm(-a_stand)*(dnorm(a_stand)*(1-a_stand^2))+a_stand*dnorm(a_stand)^2))/(pnorm(-a_stand)^2))
      H2 <- sum(((pnorm(v_stand[cens])-pnorm(a_stand[cens]))*(((dnorm(v_stand[cens])*(v_stand[cens]^3-v_stand[cens]))-((dnorm(a_stand[cens])*(a_stand[cens]^3-a_stand[cens])))))+((dnorm(v_stand[cens])*v_stand[cens])-(dnorm(a_stand[cens])*a_stand[cens]))^2)/((pnorm(v_stand[cens])-pnorm(a_stand[cens]))^2))
      H3 <- 2*sum((y_stand[uncens]^2))
      H[p,p] <- H1 - H2 - H3

      #completing the matrix using the symmetric property
      H[lower.tri(H)]<- H[upper.tri(H)]
    }
  }else if(is.null(v)==TRUE){
    a_stand <- (a-X%*%theta[1:(p-1)])/exp(log_sigma) #standardized value with respect to truncated value
    y_stand <- (y-X%*%theta[1:(p-1)])/exp(log_sigma)

    if(a == -Inf){
      #calculating the hessian for each beta parameter
      for(j in 1:(p-1)){
        for(k in j:(p-1)){
          H1 <- 0
          H3 <- (X[,j]%*%X[,k])*exp(-2*log_sigma)
          H[j, k] <- H1 - H3
        }
      }

      #calculating the hessian for beta with respect to log_sigma
      for(j in 1:(p-1)){
        H1 <- 0
        H3 <- -2*exp(-log_sigma)*(X[,j]%*%y_stand)
        H[j, p] <- H1 + H3
      }
      #this completes the upper triangle of the hessian matrix

      #calculating the hessian with respect to log_sigma

      #calculating the second derivative with respect to log_sigma squared
      H1 <- 0
      H3 <- 2*sum((y_stand^2))
      H[p,p] <- H1 - H3

      #completing the matrix using the symmetric property
      H[lower.tri(H)]<- H[upper.tri(H)]
    }else{
      #calculating the hessian for each beta parameter
      for(j in 1:(p-1)){
        for(k in j:(p-1)){
          H1 <- sum((-X[,j]*X[,k]*exp(-2*log_sigma))*((pnorm(-a_stand)*dnorm(a_stand)*a_stand-dnorm(a_stand)^2)/pnorm(-a_stand)^2))
          H3 <- (X[,j]%*%X[,k])*exp(-2*log_sigma)
          H[j, k] <- H1 - H3
        }
      }

      #calculating the hessian for beta with respect to log_sigma
      for(j in 1:(p-1)){
        H1 <- sum((X[,j]*exp(-log_sigma))*((pnorm(-a_stand)*(dnorm(a_stand)*(1-a_stand^2))+a_stand*dnorm(a_stand)^2)/pnorm(-a_stand)^2))
        H3 <- -2*exp(-log_sigma)*(X[,j]%*%y_stand)
        H[j, p] <- H1 + H3
      }
      #this completes the upper triangle of the hessian matrix

      #calculating the hessian with respect to log_sigma

      #calculating the second derivative with respect to log_sigma squared
      H1 <- sum((a_stand*(pnorm(-a_stand)*(dnorm(a_stand)*(1-a_stand^2))+a_stand*dnorm(a_stand)^2))/(pnorm(-a_stand)^2))
      H3 <- 2*sum((y_stand^2))
      H[p,p] <- H1 - H3

      #completing the matrix using the symmetric property
      H[lower.tri(H)]<- H[upper.tri(H)]
    }
  }
  return(H)
}
