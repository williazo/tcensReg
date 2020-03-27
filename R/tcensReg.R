#'  Regression Method for Truncated Normal Distribution with Censored Data
#'
#'  This function is used to find estimates from a linear equation
#'  assuming that the underlying distribution is truncated normal
#'  and the data has subsequently been censored data.
#'
#' @param formula Object of class \code{formula} which symbolically describes the model to be fit
#' @param a Numeric scalar indicating the truncation value. Initial value is -Inf indicating no truncation
#' @param v Numeric scalar indicating the censoring value. Initially set to NULL indicating no censoring
#' @param data Data.frame that contains the outcome and corresponding covariates. If none is provided then assumes objects are in user's environment.
#' @param method Character value indicating which optimization routine to perform.
#' Choices include \code{Newton}, \code{BFGS}, and \code{CG}. See details for explanation on each method.
#' @param ... Additional arguments from such as \code{max_iter}, \code{step_max}, or \code{epsilon}.
#' See details for how to define these additional arguments.
#'
#' @details
#'  This estimation procedure returns maximum likelihood estimates under the presence
#'  of \emph{left} truncation and/or \emph{left} censoring. It builds upon currently available methods
#'  for limited dependent variables by relaxing the assumption of a latent normal distribution
#'  and instead allows the this underlying distribution to potentially be a
#'  latent \strong{truncated} normal distribution.
#'
#'  To indicate left censoring the user should specify the parameter \code{v}, and
#'  to indicate left truncation specify the parameter \code{a}. If specifying
#'  both left censoring and left truncation note that there is an implicit restriction that
#'  \eqn{a<\nu}.
#'
#'  Below is a brief description of the types of distributions that can be fit
#'  along with the assumed data generating process for the observed outcome, \eqn{Y}.
#'
#'  \subsection{Latent Distributions}{
#'
#'  The tcensReg function allows user to specify one of four combinations of
#'  distributional assumptions with or without censoring. These are listed below
#'  along with the necessary arguments needed to fit this model.
#'    \subsection{Truncated Normal with Censoring}{
#'
#'    This is the main model that this package is designed to fit and introduced
#'    in \insertCite{williams2019modeling}{tcensReg}.It assumes
#'    \deqn{Y_{i}^{*}\sim TN(\mu, \sigma^{2}, a)}
#'    where TN indicates a left truncated normal random variable, truncated at
#'    the value \eqn{a}.
#'
#'    This underlying truncated normal random variable is then \emph{left} censored at the value
#'    \eqn{\nu} to create the censored observations \eqn{Y} such that
#'    \deqn{Y_{i}=\nu 1_{\{Y_{i}^{*}\le\nu\}}
#'    + Y_{i}^{*}
#'    1_{\{Y_{i}^{*}>\nu}\}}
#'
#'    Required Arguments:
#'    \itemize{
#'      \item \code{a}: left truncation value
#'      \item \code{v}: left censoring value
#'      }
#'    }
#'    \subsection{Normal with Censoring}{
#'
#'    This model is commonly referred to as the Tobit model
#'    \insertCite{tobin1958estimation}{tcensReg}. This model assumes that the
#'    data is generated from a latent normal random variable \eqn{Y_{i}^{*}}, i.e.,
#'    \deqn{Y_{i}^{*}\sim N(\mu, \sigma^{2})}
#'
#'    This underlying normal random variable is then \emph{left} censored at the value
#'    \eqn{\nu} to create the censored observations \eqn{Y} such that
#'    \deqn{Y_{i}=\nu 1_{\{Y_{i}^{*}\le\nu\}}
#'    + Y_{i}^{*}
#'    1_{\{Y_{i}^{*}>\nu}\}}
#'
#'    Required Arguments:
#'    \itemize{
#'     \item \code{v}: left censoring value
#'    }
#'
#'    This procedure can also be fit using the \code{\link[censReg]{censReg}} package by
#'    \insertCite{henningsen2010estimating;textual}{tcensReg}.
#'    }
#'    \subsection{Truncated Normal}{
#'
#'      This model assumes that there is no censored observations, but that the
#'      data are left truncated as originally described by
#'      \insertCite{hald1949maximum;textual}{tcensReg}.
#'
#'      Therefore, we assume that the observed values follow
#'      \deqn{Y_{i}^\sim TN(\mu, \sigma^{2}, a)}
#'      where TN indicates a left truncated normal random variable, truncated at
#'      the value \eqn{a}.
#'
#'      Required Arguments:
#'      \itemize{
#'       \item \code{a}: left truncation value
#'      }
#'
#'      This procedure can also be fit using the \code{\link[truncreg]{truncreg}} package by
#'      \insertCite{croissant2018truncreg;textual}{tcensReg}.
#'     }
#'     \subsection{Normal}{
#'
#'     This model assumes that there is no left censoring and no left truncation.
#'
#'     Maximum likelihood estimates are returned based on the assumption that the
#'     random variable follows the distribution
#'     \deqn{Y_{i}^\sim N(\mu, \sigma^{2})}
#'
#'     Required Arguments: None
#'
#'     This procedure can also be fit using the command \code{\link{lm}} in base R.
#'     }
#'  }
#'
#'  \subsection{Optimization Methods}{
#'
#'  Currently available optimization routines include conjugate gradient
#'  (\code{CG}), Newton-Raphson (\code{Newton}), and BFGS (\code{BFGS}).
#'  The default method is set as the conjugate gradient. Both the of the
#'  conjugate gradient and BFGS methods are implemented via the
#'  general-purpose optimization routine \code{\link{optim}}. These two methods
#'  use only the respective likelihood and gradient functions.
#'  The Newton-Raphson method uses the likelihood, gradient, and Hessian
#'  functions along with line search to achieve the maximum likelihood.
#'  }
#'
#'  \subsection{Additional Arguments}{
#'
#'  There are additional arguments that the user may provide for controlling the
#'  optimization routine.
#'  \itemize{
#'    \item{max_iter: Maximum number of iterations for optimization routine. Default is 100}
#'    \item{step_max: Maximum number of steps when performing line search. Default is 10}
#'    \item{epsilon: Numeric value used to define algorithm stops, i.e., when evaluated gradient is less than epsilon. Default is 0.001}
#'    \item{tol_val: Tolerance value used to stop the algorithm if the (n+1) and (n) log likelihood is within the tolerance limit.}
#'    \item{theta_init: Numeric vector specifying the initial values to use for the estimated parameters \eqn{\beta} and \eqn{\log(\sigma)}}
#'  }
#'  }
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
#' trunc_cens_mod <- tcensReg(y ~ 1, v = 0.25, a = 0)
#' summary(trunc_cens_mod)
#'
#'
#' @return Returns a list of final estimate of theta, total number of iterations performed, initial log-likelihood,
#' final log-likelihood, estimated variance covariance matrix, information criterion, model design matrix,
#' call, list of total observations and censored observations, and latent distributional assumption.
#'
#' @export
#' @references
#'     \insertAllCited{}
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

    #observed censoring
    n_total <- length(y)
    n_cens <- sum(y==v)

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

    #additional information to save in order to use S3 methods
    results$model_matrix <- X
    results$call <- match.call()
    results$n_count <- list(n=n_total, n_cens=n_cens)
    latent_assumption <- ifelse(is.null(v) & a == -Inf, "Normal",
                                ifelse(is.null(v) & a != -Inf, "Truncated Normal",
                                       ifelse(!is.null(v) & a == -Inf, "Normal with Censoring",
                                              ifelse(!is.null(v) & a != -Inf, "Truncated Normal with Censoring"))))
    results$latent_assumption <- latent_assumption
    #saving object as S3
    class(results) <- "tcensReg"
    return(results)
}
