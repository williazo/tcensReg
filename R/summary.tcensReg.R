#'  Method for Returning Summary Statistics for tcensReg S3 Objects
#'
#' @method summary tcensReg
#'
#' @param object Object of class \code{tcensReg}
#' @param logsigma Logical indicator of whether to return standard deviation as logsigma or sigma. Default is FALSE returning sigma
#' @param digits Numeric value for number of digits to display for all numeric results
#' @param ... Additional optional parameters
#'
#' @return Relevant results of tcensReg fitted object
#' @export
#' @keywords internal
summary.tcensReg <- function(object, logsigma=FALSE, digits=4, ...){
    print_ll <- round(object$final_ll, digits)
    r2_obj <- pseudo_r2(object, type="m")
    r2 <- round(r2_obj$r2, digits)
    est <- coef(object, logsigma, digits)
    #need to right a seperate var.cov.tcensReg function to get standard errors depending on logsigma value
    se <- sqrt(diag(vcov(object, logsigma)))
    t_val <- est/se
    coef_tbl <- cbind(est, se, t_val)
    colnames(coef_tbl) <- c("Estimate", "Std. Error", "t value")
    cat("\n")
    cat("Call:\n")
    cat(paste0(deparse(object$call), sep = "\n", collapse = "\n" ))
    cat("\n")
    cat("Assumed Distribution:\n")
    cat(object$latent_assumption)
    cat("\n\n")
    cat("Count of Observations:\n")
    cnt_tbl <- unlist(object$n_count)
    names(cnt_tbl) <- c("Total Observations", "Censored Observations")
    print(cnt_tbl)
    cat("\n\n")
    cat("Coefficients:\n")
    print(round(coef_tbl, digits))
    cat("\n")
    cat("Log Likelihood: ")
    cat(as.character(print_ll))
    cat("\n")
    cat("Information Criterion: ")
    cat(paste0("AIC=", round(object$info_criteria[1], digits), " BIC=", round(object$info_criteria[2], digits)))
    cat("\n")
    cat("Optimization Method: ")
    cat(object$method)
    cat("\n")
    cat("Psuedo R2: ")
    cat(paste0(r2, " method - ", r2_obj$type))
    cat("\n\n")
}
