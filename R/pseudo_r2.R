#' @title Pseudo R2 for tcensReg Objects
#'
#' @description Implementation of various methods for calculating pseudo R2 values
#' popular with censored observations
#'
#' @param obj Object of class \code{tcensReg}
#' @param type Character value indicating the type of pseudo R2 to calculate.
#' Currently only mckelvey_zavoina is available
#'
#' @details
#' When comparing goodness of fit between methods for censored observations
#' pseudo \eqn{R^2} is often the preferred metric
#' \insertCite{veall1996pseudo}{tcensReg}. While there are many different types
#' of pseudo \eqn{R^2} measures available this function implements those that
#' are particularly relevant for censored observations. Below is a description.
#' of the currently available methods within this function.
#'
#' \subsection{McKelvey-Zavoina}{
#'
#'   This measure of pseudo \eqn{R^2} is from
#'   \insertCite{mckelvey1975statistical;textual}{tcensReg} and is the optimal
#'   metric suggested for limited dependent variables from
#'   \insertCite{veall1996pseudo;textual}{tcensReg}. The formula is shown below
#'     \deqn{R^{2}_{mz}=
#'       \frac{\sum_{i=1}^{N}(\hat{y_{i}}-\bar{\hat{y_{i}}})^{2}}
#'       {\sum_{i=1}^{N}(\hat{y_{i}}-\bar{\hat{y_{i}}})^{2}+N\hat{\sigma}}}
#' }
#'
#' @examples
#' #truncated normal underlying data
#' y_star <- rtnorm(n = 1000, mu = 0.5, sd = 1, a = 0)
#'
#' #apply censoring
#' y <- ifelse(y_star <= 0.25, 0.25, y_star)
#'
#' #find MLE estimates
#' mod_result <- tcensReg(y ~ 1, v = 0.25, a = 0)
#'
#' pseudo_r2(mod_result, type="m")
#'
#' @return List with numeric value representing the pseudo \eqn{R^2} and type of pseudo \eqn{R^2} calculated
#' @importFrom Rdpack reprompt
#'
#' @references
#'     \insertAllCited{}
#'
#' @export
pseudo_r2 <- function(obj, type=c("mckelvey_zavoina")){
    type <-  match.arg(type)
    if(class(obj)!="tcensReg"){
        stop("obj must be of class tcensReg")
    }
    if(type=="mckelvey_zavoina") {
        obj_coef <- coef(obj, logsigma=FALSE)
        sigma_hat <- unname(obj_coef[grep("sigma", names(obj_coef))])
        y_hat <- predict(obj)
        y_hatbar <- mean(y_hat)
        explained_var <- norm(y_hat - y_hatbar, "2")^2
        N <- length(y_hat)
        p_r2 <- explained_var / (explained_var + (N*sigma_hat))
    }
    pseudo_result <- list(r2=p_r2, type=type)
    return(pseudo_result)
}
