#' @title Simulate Random Left-Truncated Normal Distribution
#'
#' @description This function is used to generate random samples from left-truncated normal distribution with specified mean and variance. Sampling is performed by
#' first drawing from a random uniform distribution to generate c.d.f. probabilities and then the inverse density function is applied to generate observations.
#'
#' @param n Numeric scaler representing the number of observations. Must be greater than or equal to 1.
#' @param mu Numeric vector of means
#' @param sd Numeric vector of standard deviations.
#' @param a Numeric vector indicating the left-truncation value.
#'
#' @details Note that if the mean \code{mu} is specified as a vector then the standard deviation \code{sigma} must have the same length.
#'
#' @examples
#' #zero truncated normal data with mean 0.5 and standard deviation 1
#' y_star <- rtnorm(n = 100, mu = 0.5, sd = 1, a = 0)
#'
#' @importFrom stats runif qnorm pnorm
#'
#' @return Returns a vector of samples drawn from the specified distribition.
#'
#' @export
rtnorm <- function(n, mu, sd, a){
    probs <- stats::runif(n, 0, 1)
    if(length(mu)>1){
        if(sum(length(mu)==n, length(sd)==n) != 2) stop("Mean and standard deviation should be length 1 or length equal to n", call.=FALSE)
    } else{
        mu <- rep(mu, n)
        sd <- rep(sd, n)
    }
    if(length(a)>1) stop("Truncation is assumed to be constant across population therefore a must be length 1", call.=FALSE)
    a <- rep(a, n)
    a_stand <- (a - mu) / sd
    result <- qnorm((probs * (1 - pnorm(a_stand))) + pnorm(a_stand)) * sd + mu
    return(result)
}
