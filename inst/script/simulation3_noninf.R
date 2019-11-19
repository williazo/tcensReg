#################################
## Estimating Non-Inferiority Type I Error Rates for Difference in Means with Left Censored Data underlying Truncated Normal
## Version 1.3
## Created by:  Justin Williams
##              Alcon Intern, R&D
## Produced:    July-August 2018
#################################
#installing and loading the needed packages
list.of.packages <- c("devtools", "tictoc", "future.apply")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = T)
rm(new.packages, list.of.packages)

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

cens_diff_sim_noninf <- function(rand_seed, mu1_vec, non_inf_margin, sd_vec, n1, n2, B, tobit_val, a, alpha){
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

  #setting the intercept mean to be the mean of the first population
  beta_0 <- mu1_vec
  true_diff <- non_inf_margin
  #calculating the number of parameters in each vector
  num_mu1 <- length(mu1_vec)
  num_sd <- length(sd_vec)
  #creating the design matrix
  # X <- cbind(rep(1, n1+n2), c(rep(0, n1), rep(1, n2)))
  #the corresponding mu2 values based on the difference values
  beta_1 <- true_diff
  mu2_vec <- mu1_vec + true_diff
  beta <- cbind(beta_0 = beta_0, beta_1 = true_diff)
  # Xb <- X%*%t(beta)
  # max_trunc_prob <- pnorm(a, mean = min(mu2_vec), sd = max(sd_vec), lower.tail = TRUE)

  set.seed(rand_seed)
  ls_dt <- future.apply::future_lapply(1:B, function(num_reps){ #looping the function over the number of replicates
    lapply(sd_vec, function(s){ #applying over the number of different standard deviation values
      lapply(1:nrow(beta), function(x){ #each row of beta matrix represents a unique mu1, mu2 combination
        y_star_1 <- rtnorm(n = n1, mean = unname(beta[x, 1]), sd = s, a = a) #oversampling from a normal distribution
        y_star_2 <- rtnorm(n = n2, mean = sum(beta[x, ]), sd = s, a = a)
        y_star <- c(y_star_1, y_star_2)
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
  #dimension three is the number of different mu combinations. This is equal to length(mu1_vec)
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
  percent_cens <- data.frame(cbind(t(matrix(prop_cens, nrow = 2, ncol = num_mu1 * num_sd, dimnames = list(paste0(c("Group_1", "Group_2"), "_pctcens"), NULL))),
                                   beta_0 = mu1_vec,
                                   beta_1 = true_diff,
                                   sd = rep(sd_vec, each = num_mu1)))
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
    # comp_trunc_mod <- tcensReg(y_star ~ group, data = df, a = a)
    comp_trunc_mod <- truncreg::truncreg(y_star ~ group, data = df, point = 0)
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
    # comp_trunc_diff <- comp_trunc_mod$theta[2]
    comp_trunc_diff <- coef(comp_trunc_mod)[2]
    diff_est <- rbind(comp_diff, comp_trunc_diff, dl_diff, dl_half_diff, tobit_diff, tcensReg_diff)

    #standard deviation estimates
    comp_sd <- summary(comp_mod)$sigma
    dl_sd <- summary(dl_mod)$sigma
    dl_half_sd <- summary(dl_half_mod)$sigma
    # comp_trunc_sd <- exp(comp_trunc_mod$theta[3])
    comp_trunc_sd <- coef(comp_trunc_mod)[3]
    sd_est <- rbind(comp_sd, comp_trunc_sd, dl_sd, dl_half_sd, tobit_sd, tcensReg_sd)

    #returning the combined mean estimate and standard deviation estimates
    #first four rows are the mean estimate
    #next four rows are the standard deviation estimate
    return(rbind(diff_est, sd_est))
  })
  #these rep arrays are 12 x num_mu1 x num_sd x B
  #dim1 represent each of the 6 methods of analysis x 2 parameters (delta and sigma)
  #dim2 represents the number of different mean combinations depending on mu1
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
    non_inferior_margin <- non_inf_margin #this value is chosen based on D.8.2.2 of ISO 11979-7 2014; and C.3.2.2 of ISO 11979-9 2006 (typically -0.15 for our example)

    sd_rep_adj <- ((n1+n2)/(n1+n2-2)) * sd_rep[,non_inf_values,,] #applying an adjustment since the mle was biased
    lb_ci <- diff_rep[,non_inf_values,,] - qnorm(1 - alpha)*(sd_rep_adj*sqrt(1/n1 + 1/n2)) #calculating the lower_bound confidence interval
    reject_avg <- apply(lb_ci, c(1, 2, 3), function(x) mean(x>non_inf_margin)) #if the lb confidence interval is above the non_inf margin then we reject
    reject_df <-data.frame(t(matrix(reject_avg, nrow = 6, ncol = length(non_inf_values) * num_sd,
                                    dimnames = list(paste0(method_names, "_nullreject"), NULL))))
    reject_results <- data.frame(percent_cens[which(percent_cens$delta==-0.15), c("mu1", "mu2", "delta", "sd")], reject_df)
    writeLines("=======================================")
    writeLines("Rejecting the Null Hyptothesis of Non-Inferiority:")
    writeLines("=======================================")
    print(round(reject_results, 4))
  }

  return(reject_results)


}

tictoc::tic()
sim_results <- cens_diff_sim_noninf(rand_seed = 1509587, mu1_vec = c(1.1, 1.0), non_inf_margin = -0.15, sd_vec = c(0.4, 0.45, 0.5),
                                    n1 = 100, n2 = 100, B = 10000, tobit_val = 0.61, a = 0, alpha = 0.05)
tictoc::toc()
# parallization is done locally using 4 cores
# total time for 10000 repitions 1197.677 sec elapsed (~aprox 20 mins)
