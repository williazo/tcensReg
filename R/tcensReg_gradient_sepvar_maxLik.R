#' Gradient Vector for Truncated Normal Distribution with Censoring with Linear Equation Mean for J Independent Truncated Normals with Seprate Variance
#'
#' @param theta Numeric vector numeric vector containing estimates of beta and log sigma
#' @param left_trunc Numeric scalar indicating the truncation value
#' @param v Numeric scalar indicating the censoring value
#' @param y Numeric vector with the observed truncated and censored outcomes
#' @param X Numeric design matrix
#' @param group Character vector identifying the group membership for the independent truncated normal variables. This defines the \code{J} groups.
#'
#' @importFrom stats dnorm pnorm
#'
#' @return Vector of gradient values with p-1 beta parameters and log sigma for the nth iterate
#' @keywords internal

tcensReg_gradient_sepvar_maxLik <- function(theta, y, X, group, left_trunc = -Inf, v = NULL){

    nabla <- vector(length = length(theta)) #creating an empty vector to store the gradient

    #assume that there are a total of p parameters
    p <- length(theta)
    num_groups <- length(unique(group)) #quantifying the J groups
    log_sigmas <- theta[(p-num_groups+1):p] #these are all of the separate standard deviation parameters
    lin_pred <- theta[1:(p-num_groups)] #these are the parameters for the linear predictors

    if(is.null(v) == FALSE){ #this is the censored only version

        #calculating the gradient for each beta parameter
        for(j in 1:length(lin_pred)){
            d1 <- NULL; d2 <- NULL; d3 <- NULL; nabla_jk <- vector(length = num_groups)
            for(k in 1:num_groups){
                y_k <- y[group == unique(group)[k]]
                X_k <- as.matrix(X[group == unique(group)[k],], nrow = length(group == unique(group)[k]), ncol = ncol(X)) #I want to ensure this is still a matrix in order to use the matrix multiplicaiton later
                X_kj <- as.matrix(X[group == unique(group)[k], j], nrow = length(group == unique(group)[k]), ncol = ncol(X)) #I want to ensure this is still a matrix in order to use the matrix multiplicaiton later
                n_k <- length(y_k) #number of observations in jth group

                uncens <- which(y_k>v) #identifying which values are uncensored
                cens <- which(y_k == v) #idenfitying which values are censored

                n_0k <- length(cens) #number of censored observations
                n_1k <- length(uncens) #number of uncensored observations
                a_stand <- (left_trunc-X_k%*%lin_pred)/exp(log_sigmas[k]) #standardized value with respect to truncated value
                v_stand <- (v-X_k%*%lin_pred)/exp(log_sigmas[k]) #standardized value with respect to censored value
                y_stand <- (y_k-X_k%*%lin_pred)/exp(log_sigmas[k])


                d1 <- sum((-X_kj*exp(-log_sigmas[k])*dnorm(a_stand))/pnorm(-a_stand))
                d2 <- sum((X_kj[cens]*exp(-log_sigmas[k])*(dnorm(v_stand[cens])-dnorm(a_stand[cens])))/(pnorm(v_stand[cens])-pnorm(a_stand[cens])))
                d3 <- sum(X_kj[uncens]*exp(-log_sigmas[k])*(y_stand[uncens]))
                nabla_jk[k] <- d1 - d2 + d3
            }
            nabla[j] <- sum(nabla_jk)
        }
        k <- 1 #indicator for looping over the groups
        for(j in (p-num_groups+1):p){ #looping over the sigma parameters
            y_k <- y[group == unique(group)[k]]
            X_k <- as.matrix(X[group == unique(group)[k],], nrow = length(group == unique(group)[k]), ncol = ncol(X)) #I want to ensure this is still a matrix in order to use the matrix multiplicaiton later
            n_k <- length(y_k) #number of observations in jth group

            uncens <- which(y_k>v) #identifying which values are uncensored
            cens <- which(y_k == v) #idenfitying which values are censored

            n_0k <- length(cens) #number of censored observations
            n_1k <- length(uncens) #number of uncensored observations
            a_stand <- (left_trunc-X_k%*%lin_pred)/exp(log_sigmas[k]) #standardized value with respect to truncated value
            v_stand <- (v-X_k%*%lin_pred)/exp(log_sigmas[k]) #standardized value with respect to censored value
            y_stand <- (y_k-X_k%*%lin_pred)/exp(log_sigmas[k])

            #here we need a special case when there is no truncation (i.e. a = -Inf)
            if(left_trunc == -Inf){
                d1 <- 0
                d2 <- sum((v_stand[cens]*dnorm(v_stand[cens]))/(pnorm(v_stand[cens])))
                d3 <- n_1k
                d4 <- sum(y_stand[uncens]^2)
                nabla[j] <- d1 - d2 - d3 + d4
            }else{
                d1 <- sum((-a_stand*dnorm(a_stand))/pnorm(-a_stand))
                d2 <- sum((v_stand[cens]*dnorm(v_stand[cens])-a_stand[cens]*dnorm(a_stand[cens]))/(pnorm(v_stand[cens])-pnorm(a_stand[cens])))
                d3 <- n_1k
                d4 <- sum(y_stand[uncens]^2)
                nabla[j] <- d1 - d2 - d3 + d4
            }
            k <- k + 1
        }
    }else if(is.null(v)==TRUE){
        #calculating the gradient for each beta parameter
        for(j in 1:length(lin_pred)){
            d1 <- NULL; d2 <- NULL; d3 <- NULL; nabla_jk <- vector(length = num_groups)
            for(k in 1:num_groups){
                y_k <- y[group == unique(group)[k]]
                X_k <- as.matrix(X[group == unique(group)[k],], nrow = length(group == unique(group)[k]), ncol = ncol(X)) #I want to ensure this is still a matrix in order to use the matrix multiplicaiton later
                X_kj <- as.matrix(X[group == unique(group)[k], j], nrow = length(group == unique(group)[k]), ncol = ncol(X)) #I want to ensure this is still a matrix in order to use the matrix multiplicaiton later
                n_k <- length(y_k) #number of observations in jth group


                a_stand <- (left_trunc-X_k%*%lin_pred)/exp(log_sigmas[k]) #standardized value with respect to truncated value
                y_stand <- (y_k-X_k%*%lin_pred)/exp(log_sigmas[k])


                d1 <- sum((-X_kj*exp(-log_sigmas[k])*dnorm(a_stand))/pnorm(-a_stand))
                d3 <- sum(X_kj*exp(-log_sigmas[k])*(y_stand))
                nabla_jk[k] <- d1 + d3
            }
            nabla[j] <- sum(nabla_jk)
        }

        k <- 1 #indicator for looping over the groups
        for(j in (p-num_groups+1):p){ #looping over the sigma parameters
            y_k <- y[group == unique(group)[k]]
            X_k <- as.matrix(X[group == unique(group)[k],], nrow = length(group == unique(group)[k]), ncol = ncol(X)) #I want to ensure this is still a matrix in order to use the matrix multiplicaiton later
            n_k <- length(y_k) #number of observations in jth group

            a_stand <- (left_trunc-X_k%*%lin_pred)/exp(log_sigmas[k]) #standardized value with respect to truncated value
            y_stand <- (y_k-X_k%*%lin_pred)/exp(log_sigmas[k])

            #here we need a special case when there is no truncation (i.e. a = -Inf)
            if(left_trunc == -Inf){
                d1 <- 0
                d3 <- n_1k
                d4 <- sum(y_stand^2)
                nabla[j] <- d1 - d3 + d4
            }else{
                d1 <- sum((-a_stand*dnorm(a_stand))/pnorm(-a_stand))
                d3 <- n_1k
                d4 <- sum(y_stand^2)
                nabla[j] <- d1 - d3 + d4
            }
            k <- k + 1
        }
    }
    return(nabla)
}
