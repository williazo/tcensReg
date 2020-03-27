#' @title Simulate Random Left-Truncated Normal Distribution
#'
#' @description This function is used to generate random samples from left-truncated normal distribution with specified mean and variance.
#'
#' @param n Numeric scalar representing the number of observations. Must be greater than or equal to 1.
#' @param mu Mean value of the underlying normal random variable
#' @param sd Standard deviation of underlying normal random variable
#' @param a Numeric vector indicating the left-truncation value.
#'
#' @details
#' Our goal is to draw samples from the left truncated normal random variable
#' \eqn{Y_{i}^{*}}. We define this distribution as
#' \deqn{Y_{i}^{*}\sim TN(\mu, \sigma^{2}, a)}
#'
#' Sampling is performed by first drawing from a random variable \eqn{Z} with a
#' uniform distribution on the interval \eqn{[0, 1]} to
#' generate cumulative density probabilities, \eqn{p}. Then the inverse density function
#' of the truncated normal random variable is applied to generate our desired
#' truncated normal observations.
#'
#' This inverse truncated normal function is shown below:
#' \deqn{Y_{i}^{*}=\Phi^{-1}\Bigg\{p\times\bigg[1-\Phi\big(\frac{a-\mu}{\sigma}\big)\bigg]
#' + \Phi\big(\frac{a-\mu}{\sigma}\big)\Bigg\}\times\sigma+\mu,}
#' where \eqn{p} represents the probabilities sampled from the uniform
#' distribution.
#'
#' \subsection{Notes}{
#' \itemize{
#'   \item{If the mean, \code{mu}, is specified as a vector then the standard
#'   deviation, \code{sigma}, must have either:
#'      \enumerate{
#'        \item{same length as \code{mu}}
#'        \item{be a scalar, indicating that all samples have constant standard deviation}
#'       }
#'    }
#' }
#' }
#' @examples
#' #zero truncated normal data with mean 0.5 and standard deviation 1
#' y_star <- rtnorm(n = 100, mu = 0.5, sd = 1, a = 0)
#'
#' @importFrom stats runif qnorm pnorm
#'
#' @return Returns a vector of samples drawn from the left truncated normal
#' distribution equal to length n.
#'
#' @export
rtnorm <- function(n, mu, sd, a){
    probs <- stats::runif(n, 0, 1)
    if(length(mu)==1 & length(sd) > 1) stop("Mean is scalar and standard deviation must also be scalar", call.=FALSE)
    if(length(mu) > 1 ){
        if(max(sum(length(mu)==n, length(sd)==n), sum(length(mu)==n, length(sd)==1)) != 2){
            stop("Mean and standard deviation should be length 1 or length equal to n", call.=FALSE)
        }
    } else{
        mu <- rep(mu, n)
        sd <- rep(sd, n)
    }
    if(length(a) > 1) stop("Truncation is assumed to be constant across population therefore a must be length 1", call.=FALSE)
    a <- rep(a, n)
    a_stand <- (a - mu) / sd
    #inverse truncated normal cdf
    result <- qnorm((probs * (1 - pnorm(a_stand))) + pnorm(a_stand)) * sd + mu
    return(result)
}
