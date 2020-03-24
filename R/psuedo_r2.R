#'  Regression Method for Truncated Normal Distribution with Censored Data
#'#'
#' @param obj Object of class \code{formula} which symbolically describes the model to be fit
#' @param type Character value indicating the type of psuedo R^2 to calculate. Currently only mckelvey_zavoina is available
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
#' psuedo_r2(mod_result, type="m")
#'
#' @return List with numeric value representing the psuedo R^2 and type of psuedo R^2 calculcated
#'
#' @export

psuedo_r2 <- function(obj, type=c("mckelvey_zavoina")){
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
    psuedo_result <- list(r2=p_r2, type=type)
    return(psuedo_result)
}
