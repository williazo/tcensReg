#'  Method for Printing tcensReg S3 Objects
#'
#' @method print tcensReg
#'
#' @param x Object of class \code{tcensReg}
#' @param logsigma Logical indicator of whether to return standard deviation as logsigma or sigma. Default is FALSE returning sigma
#' @param digits Numeric scalar indicating the number of digits to return
#' @param ... Additional optional arguments to pass into the printing process
#'
#' @return default method for printing an S3 object of class tcensReg
#' @export
#' @keywords internal


print.tcensReg <- function(x, logsigma=FALSE, digits=4, ...){
    cat("Call:\n")
    cat(paste0(deparse(x$call), sep = "\n", collapse = "\n" ))
    cat("\n")
    cat("\n")
    cat("Coefficients:\n")
    print(coef(x, logsigma, digits))
}
