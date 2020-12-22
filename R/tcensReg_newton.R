#' @title Newton-Raphson Algorithm for Truncated Normal Distribution with Censoring with Linear Equation Mean
#'
#' @description Iteratively solve the optimization log likelihood problem using
#' Newton-Raphson algorithm with analytic gradient and Hessian values and step halving.
#'
#' @inheritParams tcensReg_llike
#' @param epsilon Numeric value used to define when the algorithm should stop when the gradient is less then epsilon. Default is 0.001
#' @param theta_init Initial values of theta provided by the user. If unspecified then calculates values from OLS regression
#' @param max_iter Maximum number of iterations for algorithm. Default is 100
#' @param step_max Maximum number of steps when performing line search. Default is 10
#' @param tol_val Tolerance value used to stop the algorithm if the (n+1) and (n) log likelihood is within the tolerance limit
#'
#' @importFrom stats coef dnorm lm model.frame model.matrix pnorm
#'
#' @return Returns a list of final estimate of theta, total number of iterations performed, initial log-likelihood,
#' final log-likelihood, and estimated variance covariance matrix.
#' @keywords internal


tcensReg_newton<-function(y,
                          X,
                          a = -Inf,
                          v = NULL,
                          xi = NULL,
                          b = Inf,
                          epsilon = 1e-4,
                          tol_val = 1e-6,
                          max_iter = 100,
                          step_max = 10,
                          theta_init = NULL) {

  #starting the iteration counter at one
  i <- 1

  #want to use different inital estimates depending on whether it is truncation only, censor only, or truncated and censored
  #if censored only, normal, or truncated only then use estimates from OLS
 if(is.null(theta_init) == TRUE & (a != -Inf & is.null(v) == FALSE)){
    #if censored and truncated then use initial estimates from censored only model
    cens_mod <- suppressWarnings(tcensReg(y ~ X - 1, v = v))
    theta_init <- unname(cens_mod$theta)
 } else{
   lm_mod <- lm(y ~ X - 1)
   theta_init <- c(unname(coef(lm_mod)), log(unname(summary(lm_mod)$sigma)))
  }


  theta <- theta_init #assigning the inital value for our first iterate
  p <- length(theta) #total number of parameters
  f_0 <- tcensReg_llike(theta, y, X, a, v) #calculating the initial log likelihood value
  tol_check <- 10 * tol_val #initially setting this value so that the tolerance condition will not be satisfied (guarantees at least one iteration)
  null_ll <- f_0 #saving this value to be returned at the end
  grad_vec <- tcensReg_gradient(theta, y, X, a, v) #calculating the initial gradient
  ihess_matrix <- solve(tcensReg_hess(theta, y, X, a, v)) #caculating the hessian matrix
  #here we loop through iterations until either the maximum iterations is reached or the value of the gradient is less then the specified espilon
  while(i <= max_iter & sum(abs(grad_vec)) > epsilon  & tol_check > tol_val){
    theta_pot <- theta - ihess_matrix %*% grad_vec #potential value for next iterate
    f_1 <- tcensReg_llike(theta_pot, y, X, a, v) #evaluating the log likelihoood at potential iterate
    step_counter <- 1 #setting step counter
    #it is possible that the gradient can overshoot the maximum and so we impose a line search method
    #if the updated log likelihood value is greater than the previous log likelihood then the potential iterate is accepted
    if(f_0 < f_1){
      theta <- theta_pot
    }else{
      #if the updated log likelihood value is less than or equal to previous then we need to decrease the amount that is acting on the iterate
      while(f_0 >= f_1 & step_counter < step_max){
        #replacing the potential value using step counter of 2^(-step counter)
        theta_pot <- theta - ((1 / 2) ^ step_counter) * ihess_matrix %*% grad_vec
        #caclculating the potential updated log_likelihood
        f_1 <- tcensReg_llike(theta_pot, y, X, a, v)
        #adding to the step_counter
        step_counter <- step_counter + 1
      }
      #stop the evaluation if the step_counter is equal to the maximum steps set.
      if(step_counter == step_max){warning("Line search error. Reached max step size.", call. = FALSE)}
    }

    theta <- theta_pot #updating theta
    tol_check <- abs(f_0 - f_1)
    f_0 <- f_1 #updating the log likelihood
    grad_vec <- tcensReg_gradient(theta, y, X, a, v) #updating the gradient vector
    ihess_matrix <- solve(tcensReg_hess(theta, y, X, a, v)) #updating the inverse_hess_matrix
    i <- i + 1 #adding the next iterate counter
  }

  #using the negative of the last hessian matrix as an estimate of the variance covariance matrix
  v_cov <- -ihess_matrix
  #display warning message if the maximum number of iterations is reached
  if(i == max_iter + 1) warning("Maximum iterations reached. Interpret results with caution.", call. = FALSE)

  row.names(theta) <- c(colnames(X),"log_sigma") #adding in variable names
  colnames(theta) <- "Estimate"
  row.names(v_cov) <- row.names(theta) #using the names for the variance covariance matrix
  colnames(v_cov) <- row.names(theta)

  #returning the final estimates of theta, number of iterations, inital/final log likelihood, and estimated variance covariance matrix
  return(list(theta = theta, iterations = i - 1, initial_ll = null_ll,
              final_ll = f_0, var_cov = v_cov, method="Newton"))
}
