#'  Method for Returning Covariance of Fixed Parameters for tcensReg S3 Objects
#'
#' @method vcov tcensReg
#'
#' @param object Object of class \code{tcensReg}
#' @param logsigma Logical indicator of whether to return standard deviation as logsigma or sigma. Default is FALSE returning sigma
#' @param ... Additional optional parameters
#'
#' @details
#' If \code{logsigma}=TRUE, the default variance covariance matrix is returned.
#' Otherwise if the parameter sigma is desired then \code{logsigma}=FALSE, and
#' we use the multivariate delta method to find these values as described below:
#'
#' Original Asymptotic Result:
#' \deqn{\theta=(\beta_{1}, \dots, \beta_{(p-1)}, \log(\sigma))^{T}\sim N_{p}(\hat{\theta}, \hat{\Sigma}_{\theta})}
#'
#' Transformed Parameters:
#' \deqn{g(\theta)=(\beta_{1}, \dots, \beta_{(p-1)}, \exp\{\log(\sigma)\})^{T}\sim N_{p}(g(\hat{\theta}), \nabla_{g(\theta)}^{T}\hat{\Sigma}_{\theta}\nabla_{g(\theta)})}
#'
#' The gradient vector is defined as
#' \deqn{\nabla_{g(\theta)}=diag(1_{p-1}, \log(\sigma)\exp\{\log(\sigma)\})_{p\times p}}
#'
#' @return Numeric matrix of variance covarince values
#' @export
#' @keywords internal
vcov.tcensReg <- function(object, logsigma=FALSE, ...){
    vcov_result <- object$var_cov
    dimnames(vcov_result) <- dimnames(object$var_cov)
    if(logsigma==FALSE){
        n_param <- length(object$theta)
        orig_theta <- coef(object, logsigma=TRUE)
        sigma_ind <- grep("log_sigma", names(orig_theta))
        log_sigma <- unname(orig_theta[sigma_ind])
        nabla_delta <- diag(c(rep(1, n_param-length(sigma_ind)), log_sigma*exp(log_sigma)))
        vcov_result  <- t(nabla_delta) %*% vcov_result %*% nabla_delta
        dimnames(vcov_result) <- lapply(dimnames(object$var_cov), function(x) gsub("log_", "", x))
    }
    return(vcov_result)
}
