#################################
## Difference in Means with Left Censored Data underlying Truncated Normal
## Version 1.5
## Created by:  Justin Williams
##              Alcon Intern, R&D
## Produced:    July-August 2018
#################################
#installing and loading the needed packages
.list.of.packages <- c("msm", "devtools", "tictoc", "future.apply")
.new.packages <- .list.of.packages[!(.list.of.packages %in% installed.packages()[,"Package"])]
if(length(.new.packages)) install.packages(.new.packages)
sapply(.list.of.packages, require, character.only = T)

#enabling parallel processing
future::plan(multiprocess)

#installing the package for estimating truncated with censoring from GitHub page
devtools::install_github("williazo/tcensReg")
library(tcensReg) #my own package that is used to estimate censored only, truncated only, and truncated with censoring parameters

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

# simulating data from two left censored populations

cens_diff_sim <- function(rand_seed, mu1_vec, true_diff, sd_vec, n1, n2, B, tobit_val, a, alpha){
  # rand_seed : scalar numeric value used as the argument in the random seed generator. used for reproducability of simulation results
  # mu1_vec   : vector of opbservations for the mean of Population 1
  # true_diff : vector of difference values which defines the mean difference between Population 1 and Population 2
  # sd_vec    : vector of standard deviation values
  # n1        : scalar numeric value defining the number of observations in Sample 1
  # n2        : scalar numeric value defining the number of observations in Sample 2
  # B         : scalar numeric value indicating the number of replications to use
  # tobit_val : scalar numeric value indicating the tobit threshold value where censoring is occuring
  # a         : scalar numeric value indicating the truncation value
  # alpha     : scalar numeric value used to set the Type I probability for the non-inferiority test. Error rate should be alpha/2

  set.seed(rand_seed)
  #setting the intercept mean to be the mean of the first population
  beta_0 <- mu1_vec
  #calculating the number of parameters in each vector
  num_mu1 <- length(mu1_vec)
  num_sd <- length(sd_vec)
  num_diff <- length(true_diff)
  #creating the design matrix
  X <- cbind(rep(1, n1+n2), c(rep(0, n1), rep(1, n2)))
  #the corresponding mu2 values based on the difference values
  beta_1 <- true_diff
  beta <- cbind(beta_0 = rep(beta_0, each = num_diff), beta_1 = rep(true_diff, num_mu1))
  Xb <- X%*%t(beta)
  ls_dt <- future.apply::future_lapply(1:B, function(B){ #looping the function over the number of replicates
    lapply(sd_vec, function(s){ #applying over the number of different standard deviation values
      lapply(1:ncol(Xb), function(x){ #each column of Xb represents a unique mu1, mu2 combination
        y_star = msm::rtnorm(n = n1 + n2, mean = Xb[, x], sd = s, lower = a)
        y_dl <- cens_method(y_star, method = "DL", tobit_val)
        y_dl_half <- cens_method(y_star, method = "DL_half", tobit_val)
        y_tobit <- cens_method(y_star, method = "Tobit", tobit_val)
        cens_ind <- ifelse(y_tobit == tobit_val, 1, 0)
        group <- c(rep(0, n1), rep(1, n2))
        data.frame(y_star, y_dl, y_dl_half, y_tobit, cens_ind, group)
      })})})
  dt <- array(unlist(ls_dt), dim = c(n1+n2, 6, nrow(beta), num_sd, B),
              dimnames = list(NULL, c("y_star","y_dl", "y_dl_half", "y_tobit", "cens_ind", "group"), NULL, NULL, NULL))
  #dt is a large array containing all of the different vector values
  #dimensions are (n1+n2) x dt_vars x num_mu_combo x num_sd x num_reps
  #dimension one is the number of observations drawn in the total sample which is equal to n1+n2
  #dimension two is the number of variables in the data.frame which includes the outcome variable y for each method, censoring indicator, and group variable
  #dimension three is the number of different mu combinations. This is equal to length(mu1_vec) * length(true_diff)
  #dimension four is the number of standard deviation values
  #dimension five is the total number of replicates

  #### Calculate Censoring Percentage #####
  #calculate the percent censoring within each group over the number of mu combinations, sd combinations, and number of replicates
  prop_cens_rep <- apply(dt, c(3, 4, 5), function(X){
    df <- data.frame(X)
    cens_tbl <- round(prop.table(table(df$cens_ind, df$group), 2), 4)
    return(cens_tbl)
  })

  #averaging the censoring results over the number of replicates
  prop_cens <- apply(prop_cens_rep[c(2,4),,,], c(1,2,3), mean)
  row.names(prop_cens) <- c("Group_1", "Group_2")
  percent_cens <- data.frame(cbind(t(matrix(prop_cens, nrow = 2, ncol = num_mu1 * num_diff * num_sd, dimnames = list(paste0(c("Group_1", "Group_2"), "_pctcens"), NULL))),
                                   beta_0 = rep(mu1_vec, each = num_diff),
                                   beta_1 = rep(true_diff, times = num_mu1),
                                   sd = rep(sd_vec, each = num_diff*num_mu1)))
  percent_cens$mu1 <- percent_cens$beta_0
  percent_cens$mu2 <- percent_cens$beta_0 + percent_cens$beta_1
  percent_cens$delta <- percent_cens$beta_1

  #calculating the censoring percentage for each scenario
  writeLines("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
  writeLines("Censoring Percentage")
  writeLines("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
  print(round(percent_cens, 4))
  writeLines("")
  writeLines("")

  ##### Estimating the Parameters for each of the Six Methods
  method_names <- c("Uncens_NT", "GS", "DL", "DL_half", "Tobit", "tcensReg")
  param_rep <- future.apply::future_apply(dt, c(3, 4, 5), function(dt_X){
    df <- data.frame(dt_X)
    #fitting the models
    comp_mod <- lm(y_star ~ group, data = df)
    dl_mod <- lm(y_dl ~ group, data = df)
    dl_half_mod <- lm(y_dl_half ~ group, data = df)
    comp_trunc_mod <- tcensReg(y_star ~ group, data = df, a = a)
    if(sum(df$cens_ind == 1) == 0){ #checking  to see if there are zero censored obsrvations
      tobit_mod <- lm(y_tobit ~ group, data = df)
      tobit_diff <- coef(tobit_mod)["group"]
      tobit_sd <- summary(tobit_mod)$sigma

      #if there are no censored observations then a truncated regression should be run
      tcensReg_mod <- tcensReg(y_tobit ~ group, data = df, a = a)
      tcensReg_diff <- tcensReg_mod$theta[2]
      tcensReg_sd <- exp(tcensReg_mod$theta[3])

    } else{
      #note that for censReg it returns an estimate of the log_sd since it calculates the score equations with respect to this parameter
      tobit_mod <- tcensReg(y_tobit ~ group, data = df, v = tobit_val)
      tobit_diff <- tobit_mod$theta[2]
      tobit_sd <- exp(tobit_mod$theta[3])

      tcensReg_mod <- tcensReg(y_tobit ~ group, a = a, v = tobit_val, data = df)
      tcensReg_diff <- tcensReg_mod$theta[2]
      tcensReg_sd <- exp(tcensReg_mod$theta[3])
    }


    #difference estimates
    comp_diff <- coef(comp_mod)[2]
    dl_diff <- coef(dl_mod)[2]
    dl_half_diff <- coef(dl_half_mod)[2]
    comp_trunc_diff <- comp_trunc_mod$theta[2]
    diff_est <- rbind(comp_diff, comp_trunc_diff, dl_diff, dl_half_diff, tobit_diff, tcensReg_diff)

    #standard deviation estimates
    comp_sd <- summary(comp_mod)$sigma
    dl_sd <- summary(dl_mod)$sigma
    dl_half_sd <- summary(dl_half_mod)$sigma
    comp_trunc_sd <- exp(comp_trunc_mod$theta[3])
    sd_est <- rbind(comp_sd, comp_trunc_sd, dl_sd, dl_half_sd, tobit_sd, tcensReg_sd)

    #returning the combined mean estimate and standard deviation estimates
    #first four rows are the mean estimate
    #next four rows are the standard deviation estimate
    return(rbind(diff_est, sd_est))
  })
  #these rep arrays are 12 x (num_mu1 x num_diff) x num_sd x B
  #dim1 represent each of the 6 methods of analysis x 2 parameters (delta and sigma)
  #dim2 represents the number of different mean combinations depending on mu1 and difference vector
  #dim3 represents the number of different standard deviation parameter settings
  #dim4 represents the number of replications
  diff_rep <- param_rep[1:6,,,]
  row.names(diff_rep) <- method_names
  sd_rep <- param_rep[7:12,,,]
  row.names(sd_rep) <- method_names

  #####################################
  #####Non-Inferiority Testing#########
  #####################################
  if(sum(true_diff==-0.15)>=1){ #if the true difference is set to -0.15 then we want to run non-inferiority testing
    non_inf_values <- which(beta[, "beta_1"]== -0.15)
    non_inferior_margin <- -0.15 #this value is chosen based on D.8.2.2 of ISO 11979-7 2014; and C.3.2.2 of ISO 11979-9 2006


    lb_ci <- diff_rep[,non_inf_values,,] - qnorm(1 - alpha)*(sd_rep[,non_inf_values,,]*sqrt(1/n1 + 1/n2)) #calculating the lower_bound confidence interval
    reject_null <- ifelse(lb_ci > non_inferior_margin, 1, 0) #reject the null hypothesis of non-inferiority if the lower bound is greater than the non-inferiority margin
    reject_avg <- apply(reject_null, c(1, 2, 3), mean) #here we are calculating the percent of times where the null hypothesis was rejected
    reject_df <-data.frame(t(matrix(reject_avg, nrow = 6, ncol = length(non_inf_values) * num_sd,
                                    dimnames = list(paste0(method_names, "_nullreject"), NULL))))
    reject_results <- data.frame(percent_cens[which(percent_cens$delta==-0.15), c("mu1", "mu2", "delta", "sd")], reject_df)
    writeLines("=======================================")
    writeLines("Rejecting the Null Hyptothesis of Non-Inferiority:")
    writeLines("=======================================")
    print(round(reject_results, 4))
  }

  ######################################
  ######### Mean Parameter #############

  diff_array <- apply(diff_rep, c(1,2,3), mean) #averaging over all of the replications
  diff_df <- data.frame(t(matrix(diff_array, nrow = 6, ncol = num_diff * num_mu1 * num_sd,
                                 dimnames = list(paste0(method_names, "_est"), NULL))))
  results_diff <- data.frame(percent_cens[, c("mu1", "mu2", "delta", "sd")], diff_df)
  writeLines("++++++++++++++++++++++++++++++++++++++")
  writeLines("Estimates: Mean Difference:")
  writeLines("++++++++++++++++++++++++++++++++++++++")
  print(round(results_diff, 4))

  #sum squared bias calculation
  true_array <- array(matrix(rep(rep(true_diff, times = num_mu1), 6), nrow = 6, ncol = num_mu1 * num_diff, byrow = TRUE), dim = dim(diff_rep))
  ssb_rep <- (diff_rep - true_array)^2
  ssb_array <- apply(ssb_rep, c(1,2,3), sum)
  ssb_df <- data.frame(t(matrix(ssb_array, nrow = 6, ncol = num_diff * num_mu1 * num_sd,
                                dimnames = list(paste0(method_names, "_ssb"), NULL))))
  results_ssb <- data.frame(percent_cens[, c("mu1", "mu2", "delta", "sd")], ssb_df)
  writeLines("++++++++++++++++++++++++++++++++++++++")
  writeLines("Sum of Squared Bias: Mean Difference:")
  writeLines("++++++++++++++++++++++++++++++++++++++")
  print(round(results_ssb, 4))

  #mean squared error calculation
  mse_array <- ssb_array/B #scaling by the number of replications
  mse_df <- data.frame(t(matrix(mse_array, nrow = 6, ncol = num_diff * num_mu1 * num_sd,
                                dimnames = list(paste0(method_names, "_mse"), NULL))))
  results_mse <- data.frame(percent_cens[, c("mu1", "mu2", "delta", "sd")], mse_df)
  writeLines("++++++++++++++++++++++++++++++++++++++")
  writeLines("Mean Squared Error: Mean Difference:")
  writeLines("++++++++++++++++++++++++++++++++++++++")
  print(round(results_mse, 4))


  #Absolute Bias
  results_absb <- data.frame(sapply(which(grepl("_est", names(results_diff))), function(col_num) abs(results_diff[, col_num] - results_diff$delta)))
  names(results_absb) <- paste0(method_names, "_absb")
  absb_print <- data.frame(percent_cens[, c("mu1", "mu2", "delta", "sd")], round(results_absb, 4))
  writeLines("***************************************")
  writeLines("Absolute Bias: Mean Difference:")
  writeLines("***************************************")
  print(absb_print)

  #Percent Bias
  results_pctb <- data.frame(sapply(which(grepl("_est", names(results_diff))), function(col_num) abs(results_diff[, col_num] - results_diff$delta)/abs(results_diff$delta)))
  names(results_pctb) <- paste0(method_names, "_pctb")
  #cannot estimate percent bias when true difference is equal to zero
  results_pctb[results_pctb==Inf] <- NA
  pctb_print <- data.frame(percent_cens[, c("mu1", "mu2", "delta", "sd")], round(results_pctb, 4))
  writeLines("***************************************")
  writeLines("Percent Bias: Mean Difference:")
  writeLines("***************************************")
  print(pctb_print)

  ######################################
  ######### Standard Deviation Parameter #############
  sd_array <- apply(sd_rep, c(1,2,3), mean)
  sd_df <- data.frame(t(matrix(sd_array, nrow = 6, ncol = num_diff * num_mu1 * num_sd,
                               dimnames = list(paste0(method_names, "_est"), NULL))))
  results_sd <- data.frame(percent_cens[, c("mu1", "mu2", "delta", "sd")], sd_df)
  writeLines("++++++++++++++++++++++++++++++++++++++")
  writeLines("Estimates: Standard Deviation")
  writeLines("++++++++++++++++++++++++++++++++++++++")
  print(round(results_sd, 4))

  #sum squared bias calculation
  true_array_sd <- array(rep(rep(sd_vec, each = 6 * num_mu1 * num_diff), B), dim = dim(diff_rep))
  ssb_rep_sd <- (sd_rep - true_array_sd)^2
  ssb_array_sd <- apply(ssb_rep_sd, c(1,2,3), sum)
  ssb_df_sd <- data.frame(t(matrix(ssb_array_sd, nrow = 6, ncol = num_diff * num_mu1 * num_sd,
                                   dimnames = list(paste0(method_names, "_ssb"), NULL))))
  results_ssb_sd <- data.frame(percent_cens[, c("mu1", "mu2", "delta", "sd")], ssb_df_sd)
  writeLines("++++++++++++++++++++++++++++++++++++++")
  writeLines("Sum of Squared Bias: Standard Deviation")
  writeLines("++++++++++++++++++++++++++++++++++++++")
  print(round(results_ssb_sd, 4))

  #mean squared error calculation
  mse_array_sd <- ssb_array_sd/B #scaling by the number of replications
  mse_df_sd <- data.frame(t(matrix(mse_array_sd, nrow = 6, ncol = num_diff * num_mu1 * num_sd,
                                   dimnames = list(paste0(method_names, "_mse"), NULL))))
  results_mse_sd <- data.frame(percent_cens[, c("mu1", "mu2", "delta", "sd")], mse_df_sd)
  writeLines("++++++++++++++++++++++++++++++++++++++")
  writeLines("Mean Squared Error: Standard Deviation")
  writeLines("++++++++++++++++++++++++++++++++++++++")
  print(round(results_mse_sd, 4))


  #Absolute Bias
  results_absb_sd <- data.frame(sapply(which(grepl("_est", names(results_sd))), function(col_num) abs(results_sd[, col_num] - results_sd$sd)))
  names(results_absb_sd) <- paste0(method_names, "_absb")
  absb_print_sd <- data.frame(percent_cens[, c("mu1", "mu2", "delta", "sd")], round(results_absb_sd, 4))
  writeLines("***************************************")
  writeLines("Absolute Bias: Standard Deviation")
  writeLines("***************************************")
  print(absb_print_sd)

  #Percent Bias
  results_pctb_sd <- data.frame(sapply(which(grepl("_est", names(results_sd))), function(col_num) abs(results_sd[, col_num] - results_sd$sd)/abs(results_sd$sd)))
  names(results_pctb_sd) <- paste0(method_names, "_pctb")

  pctb_print_sd <- data.frame(percent_cens[, c("mu1", "mu2", "delta", "sd")], round(results_pctb_sd, 4))
  writeLines("***************************************")
  writeLines("Percent Bias: Standard Deviation")
  writeLines("***************************************")
  print(pctb_print_sd)


  #collecting all of the results together
  results_mean <- data.frame(cbind(percent_cens[, c("mu1", "mu2", "delta", "sd")],
                                   percent_cens[, paste0(c("Group_1", "Group_2"), "_pctcens")],
                                   diff_df, results_absb, results_pctb, ssb_df, mse_df))

  results_sd <- data.frame(cbind(percent_cens[, c("mu1", "mu2", "delta", "sd")],
                                 percent_cens[, paste0(c("Group_1", "Group_2"), "_pctcens")],
                                 sd_df, results_absb_sd, results_pctb_sd, ssb_df_sd, mse_df_sd))
  if(sum(true_diff==-0.15)>=1){
    results_noninf <- data.frame(cbind(percent_cens[, c("mu1", "mu2", "delta", "sd")],
                                       percent_cens[, paste0(c("Group_1", "Group_2"), "_pctcens")],
                                       reject_df))
    return(list(mean_diff = results_mean, sd = results_sd, non_inf = results_noninf))
  } else{
    return(list(mean_diff = results_mean, sd = results_sd))
  }
}

tictoc::tic()
sim_results <- cens_diff_sim(rand_seed = 1616, mu1_vec = c(1.1, 1.0),
                             true_diff = c(-0.3, -0.2, -0.1, 0), sd_vec = c(0.4, 0.45, 0.5),
                             n1 = 100, n2 = 100, B = 10000, tobit_val = 0.61, a = 0, alpha = 0.05)
tictoc::toc()
#For 1000 iterations this took 400.187 sec elapsed (~ aprox 7 mins)


