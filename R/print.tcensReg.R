#'  Method for Printing tcensReg S3 Objects
#''
#' @param obj Object of class \code{tcensReg}
#' @param logsigma Logical indicator of whether to return standard deviation as logsigma or sigma. Default is FALSE returning sigma
#' @param digits Numeric scalar indicating the number of digits to return
#' @param ... Additional optional arguments to pass into the printing process
#'
#' @return default method for printing an S3 object of class tcensReg
#' @export
#' @keywords internal


print.tcensReg <- function(obj, logsigma=FALSE, digits=4, ...){
    print_ll <- round(obj$final_ll, digits)
    r2_obj <- psuedo_r2(obj, type="m")
    r2 <- round(r2_obj$r2, digits)
    cat("\n")
    cat("Coefficients:\n")
    print(coef(obj, logsigma, digits))
    cat("\n")
    cat("Log Likelihood: ")
    cat(as.character(print_ll))
    cat("\n")
    cat("Information Criterion: ")
    cat(paste0("AIC=", round(obj$info_criteria[1], digits), " BIC=", round(obj$info_criteria[2], digits)))
    cat("\n")
    cat("Optimization Method: ")
    cat(obj$method)
    cat("\n")
    cat("Psuedo R2: ")
    cat(paste0(r2, " method - ", r2_obj$type))
    cat("\n\n")
}
