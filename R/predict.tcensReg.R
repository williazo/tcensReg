#'  Method for Obtaining Fitted Predictions for tcensReg S3 Objects
#'
#' @method predict tcensReg
#'
#' @param object Object of class \code{tcensReg}
#' @param newdata New values to obtain predictions
#'
#' @importFrom stats predict
#'
#' @return numeric vector of fitted values, i.e., yhat
#' @export
#' @keywords internal
predict.tcensReg <- function(object, newdata=NULL, ...){
    design_mat <- object$model_matrix
    fixed_coef <- object$theta[-grep("sigma", row.names(object$theta))]
    if(is.null(newdata)){
        y_hat <- design_mat %*% fixed_coef
    } else {
        stop("newdata prediction for tcensReg class is under construction", call.=FALSE)
    }

    return(y_hat)
}
