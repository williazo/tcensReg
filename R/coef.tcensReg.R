#'  Method for Extracting Coefficient Estimates for tcensReg S3 Objects
#'
#' @param obj Object of class \code{tcensReg}
#' @param logsigma Logical indicator of whether to return standard deviation as logsigma or sigma. Default is FALSE returning sigma
#' @param digits Numeric scalar indicating the number of digits to return
#'
#' @return numeric named vector of coefficient estimates
#' @export
#' @keywords internal


coef.tcensReg <- function(obj, logsigma=FALSE, digits=4){
    theta_result <- c(obj$theta)
    names(theta_result) <- row.names(obj$theta)
    if(logsigma==FALSE){
        sigma_ind <- grep("log_sigma", names(theta_result))
        theta_result[sigma_ind] <- exp(theta_result[sigma_ind])
        names(theta_result)[sigma_ind] <- "sigma"
    }
    theta_result <- round(theta_result, digits)
    return(theta_result)
}
