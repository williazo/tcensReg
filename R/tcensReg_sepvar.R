#' Regression Method for Truncated Normal Distribution with Censoring for Indepentend Truncated Normal Groups
#'
#' This function is used to find estimates from a linear equation assuming that the data is observed from a truncated
#'  distribution with left censoring. It uses analytically derived values of the gradient vector and hessian matrix to
#'  iteratively solve for the maximum likelihood using Newton-Raphson methods with line search. This function can also
#'  be used with censored only, truncated only, or uncensored and untruncated gaussian models.
#'
#' @param formula Object of class \code{formula} which symbolically describes the model to be fit
#' @param a Numeric scalar indicating the truncation value. Initial value is -Inf indicating no truncation
#' @param v Numeric scalar indicating the censoring value. Initially set to NULL indicating no censoring
#' @param group_var Character scalar indiciating a variable in the data.frame that defines the independent groups
#' @param ... Additional arguments from \code{\link{tcensReg_newton}} such as \code{max_iter}, \code{step_max}, or \code{epsilon}.
#'
#' @importFrom stats model.frame model.matrix
#' @export
#'
#'
#' @return Returns a list of final estimate of theta, total number of iterations performed, initial log-likelihood,
#' final log-likelihood, and estimated variance covariance matrix.

tcensReg_sepvar <- function(formula, a = -Inf, v = NULL, group_var, data = sys.frame(sys.parent()), ...){

  #checks for proper specification of formula
  if(class(formula)!= "formula"){
    stop("`formula` must be a formula", call. = FALSE)
  } else if (length(formula) != 3){
    stop("`formula` must be 2-sided", call. = FALSE)
  }
  if(is.null(v) & a == -Inf){
    warning("`v` and `a` are not specified indicating no censoring and no truncation", call. = FALSE)
  }else if(is.null(v) & a!= -Inf){
    warning("`v`is not specified indicating no censoring", call. = FALSE)
  }else if(!is.null(v) & a == -Inf){
    warning("`a` is not specified indicating no truncation", call. = FALSE)
  }

  #checking for proper specification of a and v
  if(!is.null(v) & (!is.numeric(a) | !is.numeric(v))){
    stop("`a` and `v` must both be numeric", call. = FALSE)
  } else if(!is.null(v)& (length(a)!=1 | length(v)!=1)){
    stop("`a`, and `v` must both be scalars", call. = FALSE)
  }

  if(!group_var%in%names(model.frame(data))){
    stop("`group_var` must be a variable in the data.frame", call. = FALSE)
  }

  y <- model.frame(formula, data)[, 1]
  group <- data[, group_var]
  X <- model.matrix(formula, data)
  #here we have at least some explanatory variables

  #reading in the newton raphson for the truncated censored normal
  results <- tcensReg_newton_sepvar(y, X, a, v, group, ...)
  return(results)
}
