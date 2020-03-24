#' Regression Method for Truncated Normal Distribution with Censored Data
#'
#'  This function is used to find estimates from a linear equation assuming that the underlying distribution is truncated normal
#'  and the data has subsequently been censored data. It uses analytically derived values of the gradient vector and Hessian matrix to
#'  iteratively solve for the maximum likelihood using Newton-Raphson methods with step halving line search. This function can also
#'  be used with censored only (similar to \code{\link[censReg]{censReg}}), truncated only (similar to \code{\link[truncreg]{truncreg}}), or uncensored and untruncated gaussian models.
#'
#' @param formula Object of class \code{formula} which symbolically describes the model to be fit
#' @param a Numeric scalar indicating the truncation value. Initial value is -Inf indicating no truncation
#' @param v Numeric scalar indicating the censoring value. Initially set to NULL indicating no censoring
#' @param data Data.frame that contains the outcome and corresponding covariates. If none is provided then assumes objects are in user's environment.
#' @param method Character value indicating which optimization routine to perform. Choices include \code{Newton}, \code{BFGS}, and \code{CG}. See details for explanation on each method.
#' @param ... Additional arguments from \code{\link{tcensReg_newton}} such as \code{max_iter}, \code{step_max}, or \code{epsilon}.
#'
#' @details
#'  Currently available optimization routines include conjugate gradient (\code{CG}), Newton-Raphson (\code{Newton}), and BFGS (\code{BFGS}).
#'  The default method is set as the conjugate gradient. Both the of the conjugate gradient and BFGS methods are implemented via the
#'  general-purpose optimization \code{\link{optim}}. These two methods use only the respective likelihood and gradient functions.
#'  The Newton-Raphson method uses the likelihood, gradient, and Hessian functions along with line search to achieve the maximum likelihood.
#'
#' @importFrom stats model.frame model.matrix
#'
#' @examples
#' #truncated normal underlying data
#' y_star <- rtnorm(n = 1000, mu = 0.5, sd = 1, a = 0)
#'
#' #apply censoring
#' y <- ifelse(y_star <= 0.25, 0.25, y_star)
#'
#' #find MLE estimates
#' tcensReg(y ~ 1, v = 0.25, a = 0)
#'
#' @return Returns a list of final estimate of theta, total number of iterations performed, initial log-likelihood,
#' final log-likelihood, estimated variance covariance matrix, information criterion, and model design matrix.
#'
#' @export


tcensReg <- function(formula, a = -Inf, v = NULL, data = sys.frame(sys.parent()),
                     method = c("CG", "Newton", "BFGS"), ...){
    #checks for proper specification of formula
    method <- match.arg(method)
    if(class(formula) != "formula"){
        stop("`formula` must be a formula", call. = FALSE)
        } else if (length(formula) != 3){
            stop("`formula` must be 2-sided", call. = FALSE)
        }

    #outcome vector
    y <- model.frame(formula, data)[, 1]
    #design matrix
    X <- model.matrix(formula, data)

    #checking for proper specification of a and v
    if(!is.null(v) & (!is.numeric(a) | !is.numeric(v))){
        stop("`a` and `v` must both be numeric", call. = FALSE)
        } else if(!is.null(v) & (length(a) !=1 | length(v) != 1)){
            stop("`a`, and `v` must both be scalars", call. = FALSE)
            }
    # checking model specification of truncation and censoring params
    if(is.null(v) & a == -Inf){
        warning("`v` and `a` are not specified indicating no censoring and no truncation", call. = FALSE)
        } else if(is.null(v) & a != -Inf){
            if(any(y < a)){
                stop("observed values below specified truncation `a`", call.=FALSE)
            }
            warning("`v`is not specified indicating no censoring", call. = FALSE)
            } else if(!is.null(v) & a == -Inf){
                len_cens <- sum(y == v, na.rm=TRUE)
                if(len_cens == 0){
                    stop("censoring indicated but no observed censored values",
                         call.=FALSE)
                    }
                warning("`a` is not specified indicating no truncation", call. = FALSE)
                } else {
                    if(v < a) stop("censoring specified below truncation", call.=FALSE)
                    if(any(y < a)){
                        stop("observed values below specified truncation `a`", call.=FALSE)
                    }
                    len_cens <- sum(y == v, na.rm = TRUE)
                    if(len_cens == 0){
                        stop("censoring indicated but no observed censored values",
                             call.=FALSE)
                    }
                }

    #reading in the newton raphson for the truncated censored normal
    if(method == "Newton"){
        results <- tcensReg_newton(y, X, a, v, ...)
    } else if(method %in%c("BFGS", "CG")){
        results <- tcensReg_optim(y, X, a, v, method, ...)
    }

    #adding in information criterion into the results
    aic <- (2 * length(results$theta)) - (2 * results$final_ll)
    bic <- log(length(y)) * length(results$theta) - 2 * results$final_ll
    info_criteria <- c(AIC=aic, BIC=bic)
    results[["info_criteria"]] <- info_criteria
    results$model_matrix <- X

    class(results) <- "tcensReg"

    return(results)
}
