#'  Method for Obtaining Fitted Predictions for tcensReg S3 Objects
#'
#' @param obj Object of class \code{tcensReg}
#' @param newdata New values to obtain predictions
#'
#' @return numeric vector of fitted values, i.e., yhat
#'
#' @export
#' @keywords internal



predict.tcensReg <- function(obj, newdata=NULL){
    design_mat <- obj$model_matrix
    fixed_coef <- obj$theta[-grep("sigma", row.names(obj$theta))]
    if(is.null(newdata)){
        y_hat <- design_mat %*% fixed_coef
    } else {
        stop("newdata prediction for tcensReg class is under construction", call.=FALSE)
    }

    return(y_hat)
}
