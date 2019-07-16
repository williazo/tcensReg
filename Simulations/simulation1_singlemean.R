#################################
## Left Censored Data underlying Truncated Normal (Univariate estimation)
## Version 1.6
## Created by:  Justin Williams
##              Alcon Intern, R&D
## Produced:    July-August 2018
#################################
#installing and loading the needed packages
list.of.packages <- c("msm", "devtools", "tictoc", "parallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = T)
rm(new.packages, list.of.packages)

#installing the package for estimating truncated with censoring from GitHub page
devtools::install_github("williazo/tcensReg")
#my own package that is used to estimate censored only, truncated only, and truncated with censoring parameters
library(tcensReg)

#function for censoring the data based on 12.0 CPD specifications
cens_method <- function(x, method, tobit_val){
  if(method == "DL"){
    ifelse(x < 0.61, 0.61, x)
  } else if (method == "DL_half"){
    ifelse(x < 0.61, 0.305, x)
  } else if (method == "Tobit"){
    ifelse(x < 0.61, tobit_val, x)
    #want close to 0.61 since in the real data will not be able to distinguish from 0 and -1 otherwise
    #rather than using arbitrarily precise number it is easier to detect trend using 0.60
  }
}

#generating truncated normal distribution by resampling
contrast_sens_sim_tnorm <- function(rand_seed, obs, B, mu_vec, sd_vec, tobit_val, a){
  num_mu <- length(mu_vec)
  num_sd <- length(sd_vec)
  #this is issue because of the way the simulation array is created. to apply over the correct dimension each must have more than one element
  if(B <= 1){
    stop("Number of replicates must be greater than 1", call. = FALSE)
  }else if(num_mu <= 1){
    stop("Number of mu values in `mu_vec` must be greater than 1", call. = FALSE)
  }else if(num_sd <= 1){
    stop("Number of standard deviation values in `sd_vec` must be greater than 1", call. = FALSE)
  }
  set.seed(rand_seed) #setting the seed so that the simulation is reproducible

  #generating random normal data and cutting off the censored observations according to the appropriate method
  csf_dat_tnorm <- function(n, mu, sd, reps){
    parallel::mclapply(1:B, function(B){

      y_star <- msm::rtnorm(n, mu, sd, lower = a)
      y_dl   <- cens_method(y_star, "DL", tobit_val)
      y_dl_half    <- cens_method(y_star, "DL_half", tobit_val)
      y_tobit <- cens_method(y_star, "Tobit", tobit_val)
      cens_ind <- ifelse(y_tobit == tobit_val, TRUE, FALSE)
      dat <- data.frame(y_star, y_dl, y_dl_half, y_tobit, cens_ind)
      return(dat)})
  }

  #generating data for the CSF simulation parameters
  #varying the mean from 0.8 to 1.1 and the standard deviation from 0.3 to 0.5
  dt_ls <- lapply(sd_vec, function(sd_dat) lapply(mu_vec, function(mu_dat) csf_dat_tnorm(n = obs, mu = mu_dat, sd = sd_dat, reps = B)))
  #converting from the list type of data to an array
  array_test <- array(unlist(dt_ls), dim = c(obs, 5, B, num_mu, num_sd), dimnames = list(NULL, c("y_star", "y_dl", "y_dl_half", "y_tobit", "cens_ind"), NULL, NULL ,NULL))
  #this array contains all of the data for each replication and specific parameter values
  #the dimensions are obs x 5 x B x num_mu x num_sd
  #the first dimension represents the n individuals
  #the second dimension is the 5 columns of y values and censoring indicator created by the specific censoring method
  #the third dimension represents each replication
  #the fourth dimension is the mean parameters
  #the fifth dimension is the standard deviation parameters

  #calculating the sample means for each scenario across all the replicated datasets
  apply(array_test[,1,,,], c(3,4), mean)
  #calculating the sample standard deviations for each scenario across the replicated datasets
  apply(array_test[,1,,,], c(3,4), sd)
  #calculating the percent of censoring

  #calculating the censoring percentage for each scenario
  writeLines("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
  writeLines("Censoring Percentage")
  writeLines("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
  cens_p_mat <- apply(array_test[, "y_tobit",,,], c(3, 4), function(y) paste0(round(sum(y == tobit_val)/length(y) * 100), "%"))
  row.names(cens_p_mat) <- paste0("mu_", mu_vec)
  colnames(cens_p_mat) <- paste0("sd_", sd_vec)
  print(cens_p_mat)
  writeLines("")
  writeLines("")

  #### Calculating the parameter statistics
  param_rep <- apply(array_test, c(3, 4, 5), function(dt_X){
    df <- data.frame(dt_X)
    #fitting the models
    comp_mod <- lm(y_star ~ 1, data = df)
    dl_mod <- lm(y_dl ~ 1, data = df)
    dl_half_mod <- lm(y_dl_half ~ 1, data = df)
    comp_trunc_mod <- tcensReg(y_star ~ 1, data = df, a = a)
    if(sum(df$y_tobit == tobit_val) == 0){
      tobit_mod <- lm(y_tobit ~ 1, data = df)
      tobit_mu <- coef(tobit_mod)[1]
      tobit_sd <- summary(tobit_mod)$sigma

      #if there are no censored observations then a truncated regression should be run
      tcensReg_mod <- tcensReg(y_tobit ~ 1, data = df, a = a)
      tcensReg_mu <- tcensReg_mod$theta[1]
      tcensReg_sd <- exp(tcensReg_mod$theta[2])

    } else{
      #note that for censReg it returns an estimate of the log_sd since it calculates the score equations with respect to this parameter
      tobit_mod <- tcensReg(y_tobit ~ 1, data = df, v = tobit_val)
      tobit_mu <- tobit_mod$theta[1]
      tobit_sd <- exp(tobit_mod$theta[2])

      tcensReg_mod <- tcensReg(y_tobit ~ 1, a = a, v = tobit_val, data = df)
      tcensReg_mu <- tcensReg_mod$theta[1]
      tcensReg_sd <- exp(tcensReg_mod$theta[2])
    }


    #mu estimates
    comp_mu <- coef(comp_mod)
    dl_mu <- coef(dl_mod)
    dl_half_mu <- coef(dl_half_mod)
    comp_trunc_mu <- comp_trunc_mod$theta[1]
    mu_est <- rbind(comp_mu, comp_trunc_mu, dl_mu, dl_half_mu, tobit_mu, tcensReg_mu)

    #standard deviation estimates
    comp_sd <- summary(comp_mod)$sigma
    dl_sd <- summary(dl_mod)$sigma
    dl_half_sd <- summary(dl_half_mod)$sigma
    comp_trunc_sd <- exp(comp_trunc_mod$theta[2])
    sd_est <- rbind(comp_sd, comp_trunc_sd, dl_sd, dl_half_sd, tobit_sd, tcensReg_sd)

    #returning the combined mean estimate and standard deviation estimates
    #first 6 rows are the mean estimate
    #next 6 rows are the standard deviation estimate
    return(rbind(mu_est, sd_est))
  })
  #param rep is 12 x B x num_mu x num_sd
  #dimension 1 contains all of the mu and sd parameter estimates. the first 6 rows correspond to mu and final 6 correspond to sd
  #dimension 2 is the number of replications
  #dimension 3 corresponds to the different values of mu
  #dimension 4 corresponds to the different value of sd

  #6 methods consist of Complete, Complete w/ Truncation, DL, 1/2 DL, Tobit, Tobit with Truncation
  method_names <- c("Complete", "Complete_trunc", "DL", "DL_half", "Tobit", "tcensReg")

  #estimates for all of the mu parameters
  mu_rep <- param_rep[1:6,,,]
  #estimates for all of the standard deviation parameters
  sd_rep <- param_rep[7:12,,,]

  #### Mean Estimation ####
  #averaging over all of the replications
  mu_hat <- apply(mu_rep, c(1, 3, 4), mean)
  mu_df <- data.frame(t(matrix(mu_hat, nrow = 6, ncol = num_mu * num_sd,
                                 dimnames = list(paste0(method_names, "_est"), NULL))))
  mu_est <- data.frame(mu = rep(mu_vec, num_sd), sd = rep(sd_vec, each = num_mu), mu_df)
  writeLines("++++++++++++++++++++++++++++++++++++++")
  writeLines("Estimates: Mu")
  writeLines("++++++++++++++++++++++++++++++++++++++")
  print(round(mu_est, 4))

  #sum squared bias calculation
  true_array <- array(rep(mu_vec, each = 6 * B), dim = dim(mu_rep))
  ssb_rep <- (mu_rep - true_array)^2
  ssb_array <- apply(ssb_rep, c(1,3,4), sum)
  ssb_df <- data.frame(t(matrix(ssb_array, nrow = 6, ncol = num_mu * num_sd,
                                dimnames = list(paste0(method_names, "_ssb"), NULL))))
  results_ssb <- data.frame(mu_est[, c("mu", "sd")], ssb_df)
  writeLines("++++++++++++++++++++++++++++++++++++++")
  writeLines("Sum of Squared Bias: Mu")
  writeLines("++++++++++++++++++++++++++++++++++++++")
  print(round(results_ssb, 4))

  #mean squared error calculation
  mse_array <- ssb_array/B
  mse_df <- data.frame(t(matrix(mse_array, nrow = 6, ncol = num_mu * num_sd,
                                dimnames = list(paste0(method_names, "_mse"), NULL))))
  results_mse <- data.frame(mu_est[, c("mu", "sd")], mse_df)
  writeLines("++++++++++++++++++++++++++++++++++++++")
  writeLines("Mean Squared Error: Mu")
  writeLines("++++++++++++++++++++++++++++++++++++++")
  print(round(results_mse, 4))


  #Absolute Bias
  results_absb <- data.frame(sapply(which(grepl("_est", names(mu_est))), function(col_num) abs(mu_est[, col_num] - mu_est$mu)))
  names(results_absb) <- paste0(method_names, "_absb")
  absb_print <- data.frame(mu_est[, c("mu", "sd")], round(results_absb, 4))
  writeLines("***************************************")
  writeLines("Absolute Bias: Mu")
  writeLines("***************************************")
  print(absb_print)

  #Percent Bias
  results_pctb <- data.frame(sapply(which(grepl("_est", names(mu_est))), function(col_num) abs(mu_est[, col_num] - mu_est$mu)/abs(mu_est$mu)))
  names(results_pctb) <- paste0(method_names, "_pctb")
  #cannot estimate percent bias when true difference is equal to zero
  results_pctb[results_pctb==Inf] <- NA
  pctb_print <- data.frame(mu_est[, c("mu", "sd")], round(results_pctb, 4))
  writeLines("***************************************")
  writeLines("Percent Bias: Mu")
  writeLines("***************************************")
  print(pctb_print)

  #### Standard Deviation ####
  #averaging over all of the replications
  sd_hat <- apply(sd_rep, c(1, 3, 4), mean)
  sd_df <- data.frame(t(matrix(sd_hat, nrow = 6, ncol = num_mu * num_sd,
                               dimnames = list(paste0(method_names, "_est"), NULL))))
  sd_est <- data.frame(mu = rep(mu_vec, num_sd), sd = rep(sd_vec, each = num_mu), sd_df)
  writeLines("++++++++++++++++++++++++++++++++++++++")
  writeLines("Estimates: Standard Deviation")
  writeLines("++++++++++++++++++++++++++++++++++++++")
  print(round(sd_est, 4))

  #sum squared bias calculation
  true_array_sd <- array(rep(sd_vec, each = 6 * B * num_mu), dim = dim(sd_rep))
  ssb_rep_sd <- (sd_rep - true_array_sd)^2
  ssb_array_sd <- apply(ssb_rep_sd, c(1,3,4), sum)
  ssb_df_sd <- data.frame(t(matrix(ssb_array_sd, nrow = 6, ncol = num_mu * num_sd,
                                dimnames = list(paste0(method_names, "_ssb"), NULL))))
  results_ssb_sd <- data.frame(sd_est[, c("mu", "sd")], ssb_df_sd)
  writeLines("++++++++++++++++++++++++++++++++++++++")
  writeLines("Sum of Squared Bias: Standard Deviation")
  writeLines("++++++++++++++++++++++++++++++++++++++")
  print(round(results_ssb_sd, 4))

  #mean squared error calculation
  mse_array_sd <- ssb_array_sd/B #scaling by the number of replications
  mse_df_sd <- data.frame(t(matrix(mse_array_sd, nrow = 6, ncol = num_mu * num_sd,
                                   dimnames = list(paste0(method_names, "_mse"), NULL))))
  results_mse_sd <- data.frame(sd_est[, c("mu", "sd")], mse_df_sd)
  writeLines("++++++++++++++++++++++++++++++++++++++")
  writeLines("Mean Squared Error: Standard Deviation")
  writeLines("++++++++++++++++++++++++++++++++++++++")
  print(round(results_mse_sd, 4))


  #Absolute Bias
  results_absb_sd <- data.frame(sapply(which(grepl("_est", names(sd_est))), function(col_num) abs(sd_est[, col_num] - sd_est$sd)))
  names(results_absb_sd) <- paste0(method_names, "_absb")
  absb_print_sd <- data.frame(sd_est[, c("mu", "sd")], round(results_absb_sd, 4))
  writeLines("***************************************")
  writeLines("Absolute Bias: Standard Deviation")
  writeLines("***************************************")
  print(absb_print_sd)

  #Percent Bias
  results_pctb_sd <- data.frame(sapply(which(grepl("_est", names(sd_est))), function(col_num) abs(sd_est[, col_num] - sd_est$sd)/abs(sd_est$sd)))
  names(results_pctb_sd) <- paste0(method_names, "_pctb")
  #cannot estimate percent bias when true difference is equal to zero
  results_pctb_sd[results_pctb_sd==Inf] <- NA
  pctb_print_sd <- data.frame(sd_est[, c("mu", "sd")], round(results_pctb_sd, 4))
  writeLines("***************************************")
  writeLines("Percent Bias: Standard Deviation")
  writeLines("***************************************")
  print(pctb_print_sd)

  #collecting all of the results together
  results_mean <- data.frame(cbind(mu_est[, c("mu", "sd")],
                                   cens_pct = as.vector(cens_p_mat),
                                   mu_df, results_absb, results_pctb, ssb_df, mse_df))

  results_sd <- data.frame(cbind(sd_est[, c("mu", "sd")],
                                 cens_pct = as.vector(cens_p_mat),
                                 sd_df, results_absb_sd, results_pctb_sd, ssb_df_sd, mse_df_sd))
  return(list(mu = results_mean, sd = results_sd))
}

tictoc::tic()
tnorm_results <- contrast_sens_sim_tnorm(rand_seed = 123580, obs = 100, B = 1000, mu_vec = seq(0.7, 1.1, 0.1),
                                         sd_vec = c(0.4, 0.45, 0.5), tobit_val = 0.61, a = 0)
tictoc::toc()

