#'  Regression Method for Truncated Normal Distribution with Censored Data
#'
#'  This function is used to find estimates from a linear equation
#'  assuming that the underlying distribution is truncated normal
#'  and the data has subsequently been censored data.
#'
#' @param formula Object of class \code{formula} which symbolically describes
#' the model to be fit
#' @param a Numeric scalar indicating the left truncation value. Default is -Inf
#' indicating no truncation
#' @param v Numeric scalar indicating the left censoring value. Default is NULL
#' indicating no censoring
#' @param xi Numeric scalar indicating the right censoring valye. Default is
#' NULL indicating to right censoring
#' @param b Numeric scalar indicating the right truncation value. Default is Inf
#' indicating no truncation
#' @param data Data.frame that contains the outcome and corresponding
#' covariates. If none is provided then assumes objects are in user's
#' environment.
#' @param method Character value indicating which optimization routine to
#' perform. Choices include \code{Newton}, \code{BFGS}, and \code{CG}. See
#' details for explanation on each method. Default is \code{CG}.
#' @param ... Additional arguments from such as \code{max_iter},
#' \code{step_max}, or \code{epsilon}. See details for how to define these
#' additional arguments.
#'
#' @details
#'  This estimation procedure returns maximum likelihood estimates under the
#'  presence of \emph{left}/\emph{right} truncation and/or
#'  \emph{left}/\emph{right} censoring. It builds upon currently available methods
#'  for limited dependent variables by relaxing the assumption of a latent normal distribution
#'  and instead allows the this underlying distribution to potentially be a
#'  latent \strong{truncated} normal distribution.
#'
#'  To indicate left censoring the user should specify the parameter \code{v},
#'  to indicate left truncation specify the parameter \code{a}, to specify right
#'  right censoring specify \code{xi}, and right truncation the parameter
#'  \code{b}. Implicit restrictions on these values exist such that
#'  \eqn{a<\nu<\xi<b}.
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
#'    \deqn{Y_{i}^{*}\sim TN(\mu, \sigma^{2}, a, b)}
#'    where TN indicates a left truncated normal random variable, left truncated
#'    at the value \eqn{a} and right truncated at \eqn{b}.
#'
#'    This underlying truncated normal random variable is then \emph{left} censored at the value
#'    \eqn{\nu} to create the censored observations \eqn{Y} such that
#'    \deqn{Y_{i}=\nu 1_{\{Y_{i}^{*}\le\nu\}}
#'    + Y_{i}^{*}
#'    1_{\{\nu<Y_{i}^{*}<\xi}\}} + \xi 1_{\{Y_{i}^{*}\ge\xi\}}
#'
#'    Required Arguments:
#'    \itemize{
#'      \item \code{a}: left truncation value
#'      \item \code{v}: left censoring value
#'      \item \code{xi}: right censoring value
#'      \item \code{b}: right truncation value
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
#'    \eqn{\nu} and or \emph{right} censored at the value \eqn{\xi} to create
#'    the censored observations \eqn{Y} such that
#'    \deqn{Y_{i}=\nu 1_{\{Y_{i}^{*}\le\nu\}}
#'    + Y_{i}^{*}
#'    1_{\{\nu<Y_{i}^{*}<\xi}\}} + \xi 1_{\{Y_{i}^{*}\ge\xi\}}
#'
#'    Required Arguments:
#'    \itemize{
#'     \item \code{v}: left censoring value
#'     \item \code{xi}: right censoring value
#'    }
#'
#'    This procedure can also be fit using the \code{\link[censReg]{censReg}} package by
#'    \insertCite{henningsen2010estimating;textual}{tcensReg}.
#'    }
#'    \subsection{Truncated Normal}{
#'
#'      This model assumes that there is no censored observations, but that the
#'      data are left and/or right truncated as originally described by
#'      \insertCite{hald1949maximum;textual}{tcensReg}.
#'
#'      Therefore, we assume that the observed values follow
#'      \deqn{Y_{i}^\sim TN(\mu, \sigma^{2}, a, b)}
#'      where TN indicates a left truncated normal random variable, truncated at
#'      the value \eqn{a}.
#'
#'      Required Arguments:
#'      \itemize{
#'       \item \code{a}: left truncation value
#'       \item \code{b}: right truncation value
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
#' #zero left truncated normal underlying data
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
tcensReg <- function(formula,
                     a = -Inf,
                     v = NULL,
                     xi = NULL,
                     b = Inf,
                     data = sys.frame(sys.parent()),
                     method = "CG",
                     ...) {
    #checks for proper specification of formula
    method <- match.arg(method, choices = c("CG", "Newton", "BFGS"))
    if (class(formula) != "formula") {
        stop("`formula` must be a formula", call. = FALSE)
        } else if (length(formula) != 3) {
            stop("`formula` must be 2-sided", call. = FALSE)
        }

    #outcome vector
    y <- model.frame(formula, data)[, 1]
    #design matrix
    X <- model.matrix(formula, data)

    #observed censoring
    n_total <- length(y)
    #left-censoring
    n_cens_L <- sum(y == v)
    #right-censoring
    n_cens_R <- sum(y == xi)


    # checking model specification of truncation and censoring params
    cens_trunc_param_check(y = y, a = a, v = v, xi = xi, b = b)

    #reading in the newton raphson for the truncated censored normal
    if (method == "Newton") {
        results <- tcensReg_newton(y, X, a, v, xi, b, ...)
    } else if (method %in% c("BFGS", "CG")) {
        results <- tcensReg_optim(y, X, a, v, xi, b, method, ...)
    }

    #adding in information criterion into the results
    aic <- (2 * length(results$theta)) - (2 * results$final_ll)
    bic <- log(length(y)) * length(results$theta) - 2 * results$final_ll
    info_criteria <- c(AIC = aic, BIC = bic)
    results[["info_criteria"]] <- info_criteria

    #additional information to save in order to use S3 methods
    results$model_matrix <- X
    results$call <- match.call()
    results$n_count <- list(
        n = n_total,
        n_cens = sum(n_cens_L, n_cens_R),
        n_cens_L = n_cens_L,
        n_cens_R = n_cens_R
        )
    latent_assumption <- ifelse(
        is.null(v) & a == -Inf & is.null(xi) & b == Inf,
        "Normal",
        ifelse(
            is.null(v) & is.null(xi) & (a != -Inf | b != Inf),
            "Truncated Normal",
            ifelse(
                (!is.null(v) | !is.null(xi)) & a == -Inf & b == Inf,
                "Normal with Censoring",
                ifelse(
                    (!is.null(v) & (a != -Inf | b != Inf))
                    | (!is.null(xi) & (a != -Inf | b != Inf)),
                    "Truncated Normal with Censoring",
                    NA
                    )
                )
            )
        )
    results$latent_assumption <- latent_assumption
    #saving object as S3
    class(results) <- "tcensReg"
    return(results)
}

#' Checking the censored and truncation values
#'
#' @description Internal helper function to perform basic parameter checks
#'
#' @param y Observed outcome
#' @param a Left-truncation value
#' @param v Left-censoring value
#' @param xi Right-censoring value
#' @param b Right-truncation value
#'
#' @keywords internal
cens_trunc_param_check <- function(y,
                                   a,
                                   v,
                                   xi,
                                   b) {

    #identify type of censoring
    if (all(is.null(v), is.null(xi))) {
        any_cens <- FALSE
        n_cens_L <- 0
        n_cens_R <- 0
        left_cens <- FALSE
        right_cens <- FALSE
    } else {
        any_cens <- TRUE
        left_cens <- ifelse(!is.null(v), TRUE, FALSE)
        if (left_cens) {
            if (length(v) != 1 | !is.numeric(v)) {
                stop(
                    "Left censoring value must be numeric scalar",
                    call. = FALSE
                    )
            }
        }
        n_cens_L <- sum(y == v, na.rm = TRUE)
        if (n_cens_L == 0 & left_cens) {
            stop(
                "Left censoring specified but no censored observations",
                call. = FALSE
                )
        }
        right_cens <- ifelse(!is.null(xi), TRUE, FALSE)
        if (right_cens) {
            if (length(xi) != 1 | !is.numeric(xi)) {
                stop(
                    "Right censoring value must be numeric scalar",
                     call. = FALSE
                     )
            }
        }
        n_cens_R <- sum(y == xi, na.rm = TRUE)
        if (n_cens_R == 0 & right_cens) {
            stop(
                "Right censoring specified but no censored observations",
                call. = FALSE
                )
        }
        if (left_cens & right_cens) {
            if (v >= xi) {
                stop(
                    "Left censoring greater than or equal to right censoring",
                    call. = FALSE
                    )
                }
        }
    }

    #identify type of truncation
    if (all(a == -Inf, b == Inf)) {
        any_trunc <- FALSE
        left_trunc <- FALSE
        right_trunc <- FALSE
    } else {
        any_trunc <- TRUE
        if (length(a) != 1 | !is.numeric(a)) {
            stop("Left truncation value must be numeric scalar", call. = FALSE)
        }
        if (length(b) != 1 | !is.numeric(b)) {
            stop("Right truncation value must be numeric scalar", call. = FALSE)
        }
        left_trunc <- ifelse(a != -Inf, TRUE, FALSE)
        right_trunc <- ifelse(b != Inf, TRUE, FALSE)
        if (left_trunc & right_trunc & (a >= b)) {
            stop(
                "Left truncation greater than or equal to right truncation",
                call. = FALSE
                )
        }
    }

    #identify combination of censoring and truncation
    if (all(!any_cens, !any_trunc)) {
        message("No censoring and no truncation")
    } else if (left_trunc & right_trunc) {
        #truncation on both sides
        if (any(y < a, y > b)) {
            stop("Observed values outside truncation value", call. = FALSE)
        }
        if (left_cens & right_cens) {
            if ((a >= v) | (xi >= b)) {
                stop(
                    "Censoring specified after truncation values",
                    call. = FALSE
                    )
                }
            message("Left/right censoring with left/right truncation")
            } else if (left_cens & !right_cens) {
            if (a >= v) {
                stop("Censoring specified below left truncation", call. = FALSE)
            }
            message("Left censoring with left/right truncation")
            } else if (!left_cens & right_cens) {
            if (xi > b) {
                stop(
                    "Censoring specified above right truncation",
                    call. = FALSE
                    )
            }
            message("Right censoring with left/right truncation")
        } else {
            message("No censoring with left/right truncation")
        }
    } else if (left_trunc & !right_trunc) {
        if (any(y < a)) {
            stop("Observed values below left truncation value", call. = FALSE)
            }
        if (left_cens & right_cens) {
            if ((a >= v) | (a >= xi)) {
                stop("Censoring specified below/equal to left truncation", call. = FALSE)
            }
            message("Left/right censoring with left truncation")
        } else if (left_cens & !right_cens) {
            if (a >= v) {
                stop("Censoring specified below/equal to left truncation", call. = FALSE)
            }
            message("Left censoring with left truncation")
        } else if (!left_cens & right_cens) {
            if (a >= xi) {
                stop("Censoring specified below/equal to left truncation", call. = FALSE)
            }
            #unnatural condition
            warning("Right censoring with left truncation", call. = FALSE)
        } else {
            message("No censoring with left truncation")
        }
    } else if (!left_trunc & right_trunc) {
        if (any(y > b)) {
            stop("Observed values above right truncation value", call. = FALSE)
            }
        if (left_cens & right_cens) {
            if ((xi >= b) | (v >= b)) {
                stop("Censoring specified above/equal to right truncation", call. = FALSE)
            }
            message("Left/right censoring with right truncation")
        } else if (left_cens & !right_cens) {
            if (v >= b) {
                stop(
                    "Censoring specified above/equal to right truncation",
                     call. = FALSE
                    )
            }
            #unnatural condition
            warning("Left censoring with right truncation", call. = FALSE)
        } else if (!left_cens & right_cens) {
            if (xi >= b) {
                stop("Censoring specified above/equal to right truncation", call. = FALSE)
            }
            message("Right censoring with right truncation")
        } else {
            message("No censoring with right truncation")
        }
    } else {
        if (left_cens & right_cens) {
            message("Left/right censoring with no truncation")
        } else if (left_cens & !right_cens) {
            message("Left censoring with no truncation")
        } else if (!left_cens & right_cens) {
            message("Right censoring with no truncation")
        } else {
            message("No censoring with no truncation")
            }
        }
}