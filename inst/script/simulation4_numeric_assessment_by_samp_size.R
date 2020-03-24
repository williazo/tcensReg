#################################
## Left Censored Data underlying Truncated Normal (Univariate estimation)
## Numerical Assessment of MLE as Function of Sample Size
## Version 1.0
## Created by:  Justin Williams
## Produced:    March 2020
#################################
#installing and loading the needed packages
rm(list=ls())
.list.of.packages <- c("devtools", "tictoc", "future.apply", "tcensReg", "pbapply")
lapply(.list.of.packages, function(x) if(!requireNamespace(x, quietly=TRUE)) install.packages(x))
lapply(.list.of.packages, require, character.only=TRUE, quietly=TRUE)

#enabling parallel processing
future::plan(multisession)

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
contrast_sens_sim_tnorm <- function(rand_seed, obs, B, mu, sd, tobit_val, a, silent=FALSE){
    #this is issue because of the way the simulation array is created. to apply over the correct dimension each must have more than one element
    requireNamespace("future.apply")

    if(B <= 1){
        stop("Number of replicates must be greater than 1", call. = FALSE)
    }else if(length(mu) > 1){
        stop("Specify `mu` as single scalar value", call. = FALSE)
    }else if(length(sd) > 1){
        stop("Specify `sd` as single scalar value", call. = FALSE)
    }
    set.seed(rand_seed) #setting the seed so that the simulation is reproducible

    #generating random normal data and cutting off the censored observations according to the appropriate method
    csf_dat_tnorm <- function(n, mu, sd, reps){
        future.apply::future_lapply(seq_len(B), function(B){
            y_star    <- tcensReg::rtnorm(n, mu, sd, a)
            y_dl      <- cens_method(y_star, "DL", tobit_val)
            y_dl_half <- cens_method(y_star, "DL_half", tobit_val)
            y_tobit   <- cens_method(y_star, "Tobit", tobit_val)
            cens_ind  <- ifelse(y_tobit == tobit_val, TRUE, FALSE)
            dat       <- data.frame(y_star, y_dl, y_dl_half, y_tobit, cens_ind)
            return(dat)}, future.seed=rand_seed)
    }

    #generating data for the CSF simulation parameters
    dt_ls <- csf_dat_tnorm(n = obs, mu = mu, sd = sd, reps = B)
    #converting from the list type of data to an array
    array_test <- array(unlist(dt_ls), dim = c(obs, 5, B), dimnames = list(NULL, c("y_star", "y_dl", "y_dl_half", "y_tobit", "cens_ind"), NULL))
    #this array contains all of the data for each replication and specific parameter values
    #the dimensions are obs x 5 x B
    #the first dimension represents n, i.e., the number of individuals
    #the second dimension is the 5 columns of y values and censoring indicator created by the specific censoring method
    #the third dimension represents each replication

    #calculating the sample means for each scenario across all the replicated datasets using y_star
    apply(array_test[,1,], 2, mean) #these should be close to our value of mu
    #calculating the sample standard deviations for each scenario across the replicated datasets
    apply(array_test[,1,], 2, sd) #these should be close to value of sd
    #calculating the percent of censoring

    #average censoring percentage across the B replicates
    cens_p <- mean(apply(array_test[, "y_tobit",], 2, function(y) round(sum(y == tobit_val)/length(y) * 100)))

    #### Calculating the parameter statistics
    param_rep <- future.apply::future_apply(array_test, 3, function(dt_X){
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
    #param rep is 12 x B
    #dimension 1 contains all of the mu and sd parameter estimates. the first 6 rows correspond to mu and final 6 correspond to sd
    #dimension 2 is the number of replications

    #6 methods consist of Complete, Complete w/ Truncation, DL, 1/2 DL, Tobit, Tobit with Truncation
    method_names <- c("Complete", "Complete_trunc", "DL", "DL_half", "Tobit", "tcensReg")

    #estimates for all of the mu parameters
    mu_rep <- param_rep[1:6,]
    #estimates for all of the standard deviation parameters
    sd_rep <- param_rep[7:12,]

    #### Mean Estimation ####
    #averaging over all of the replications
    mu_hat <- apply(mu_rep, 1, mean)
    names(mu_hat) <- paste0(method_names)
    mu_est <- data.frame(mu = mu, est=mu_hat)

    #sum squared bias calculation
    true_array <- array(rep(mu, each = 6 * B), dim = dim(mu_rep))
    ssb_rep <- (mu_rep - true_array)^2
    ssb_array <- apply(ssb_rep, 1, sum)
    names(ssb_array) <- paste0(method_names)
    results_ssb <- data.frame(mu=mu, ssb=ssb_array)

    #mean squared error calculation
    mse_array <- ssb_array/B
    results_mse <- data.frame(mu=mu, mse=mse_array)

    #Absolute Bias
    results_absb <- abs(mu_est$est - mu_est$mu)
    names(results_absb) <- paste0(method_names)
    absb_print <- data.frame(mu=mu, absb=round(results_absb, 4))

    #Percent Bias
    results_pctb <- abs(mu_est$est - mu_est$mu)/abs(mu_est$mu)
    names(results_pctb) <- method_names
    #cannot estimate percent bias when true difference is equal to zero
    results_pctb[results_pctb==Inf] <- NA
    pctb_print <- data.frame(mu=mu, pctb=round(results_pctb, 4))

    #### Standard Deviation ####
    #averaging over all of the replications
    sd_hat <- apply(sd_rep, 1, mean)
    names(sd_hat) <- paste0(method_names)
    sd_est <- data.frame(sd = sd, est=sd_hat)

    #sum squared bias calculation
    true_array_sd <- array(rep(sd, each = 6 * B), dim = dim(sd_rep))
    ssb_rep_sd <- (sd_rep - true_array_sd)^2
    ssb_array_sd <- apply(ssb_rep_sd, 1, sum)
    names(ssb_array_sd) <- paste0(method_names)
    results_ssb_sd <- data.frame(sd=sd, ssb=ssb_array_sd)

    #mean squared error calculation
    mse_array_sd <- ssb_array_sd/B #scaling by the number of replications
    results_mse_sd <- data.frame(sd=sd, mse_array_sd)

    #Absolute Bias
    results_absb_sd <- abs(sd_est$est - sd_est$sd)
    names(results_absb_sd) <- method_names
    absb_print_sd <- data.frame(sd=sd, absb=round(results_absb_sd, 4))

    #Percent Bias
    results_pctb_sd <- abs(sd_est$est - sd_est$sd)/abs(sd_est$sd)
    names(results_pctb_sd) <- method_names
    #cannot estimate percent bias when true difference is equal to zero
    results_pctb_sd[results_pctb_sd==Inf] <- NA
    pctb_print_sd <- data.frame(sd=sd, pctb=round(results_pctb_sd, 4))

    if(silent==FALSE){
        writeLines("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        writeLines("Censoring Percentage")
        writeLines("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        print(cens_p)
        writeLines("")
        writeLines("")

        writeLines("++++++++++++++++++++++++++++++++++++++")
        writeLines("Estimates: Mu")
        writeLines("++++++++++++++++++++++++++++++++++++++")
        print(round(mu_est, 4))

        writeLines("++++++++++++++++++++++++++++++++++++++")
        writeLines("Sum of Squared Bias: Mu")
        writeLines("++++++++++++++++++++++++++++++++++++++")
        print(round(results_ssb, 4))

        writeLines("++++++++++++++++++++++++++++++++++++++")
        writeLines("Mean Squared Error: Mu")
        writeLines("++++++++++++++++++++++++++++++++++++++")
        print(round(results_mse, 4))

        writeLines("***************************************")
        writeLines("Absolute Bias: Mu")
        writeLines("***************************************")
        print(absb_print)

        writeLines("***************************************")
        writeLines("Percent Bias: Mu")
        writeLines("***************************************")
        print(pctb_print)

        writeLines("++++++++++++++++++++++++++++++++++++++")
        writeLines("Estimates: Standard Deviation")
        writeLines("++++++++++++++++++++++++++++++++++++++")
        print(round(sd_est, 4))

        writeLines("++++++++++++++++++++++++++++++++++++++")
        writeLines("Sum of Squared Bias: Standard Deviation")
        writeLines("++++++++++++++++++++++++++++++++++++++")
        print(round(results_ssb_sd, 4))

        writeLines("++++++++++++++++++++++++++++++++++++++")
        writeLines("Mean Squared Error: Standard Deviation")
        writeLines("++++++++++++++++++++++++++++++++++++++")
        print(round(results_mse_sd, 4))

        writeLines("***************************************")
        writeLines("Absolute Bias: Standard Deviation")
        writeLines("***************************************")
        print(absb_print_sd)

        writeLines("***************************************")
        writeLines("Percent Bias: Standard Deviation")
        writeLines("***************************************")
        print(pctb_print_sd)
    }

    #collecting all of the results together
    results_mean <- data.frame(cbind(mu=mu_est$mu,
                                     cens_pct = cens_p,
                                     mu_hat, results_absb, results_pctb, ssb_array, mse_array))

    results_sd <- data.frame(cbind(sd=sd_est$sd,
                                   cens_pct = cens_p,
                                   sd_hat, results_absb_sd, results_pctb_sd, ssb_array_sd, mse_array_sd))
    return(list(mu = results_mean, sd = results_sd))
}

samp_sizes <- c(100, 200, 400, 800, 1600)
tictoc::tic()
tnorm_results_list <- pbapply::pblapply(samp_sizes, function(n){
    contrast_sens_sim_tnorm(rand_seed = 123580, obs = n, B = 10000, mu = 0.9,
                            sd = 0.45, tobit_val = 0.61, a = 0, silent=TRUE)
})
tictoc::toc()
# parallization is done locally using 4 cores
# total time for 100 repitions 27.019 sec elapsed (~aprox 0.5 mins)
# total time for 10000 repitions 2279.807 sec (~approx 38 mins)

names(tnorm_results_list) <- paste0("n_", samp_sizes)

cens_diff_sim_noninf <- function(rand_seed, mu1, non_inf_margin, sd, n1, n2, B, tobit_val, a, alpha, silent=FALSE){
    # rand_seed         : scalar numeric value used as the argument in the random seed generator. used for reproducability of simulation results
    # mu1               : mean of Population 1
    # non_inf_margin    : vector of difference values which defines the mean difference between Population 1 and Population 2
    # sd                : global standard deviation
    # n1                : scalar numeric value defining the number of observations in Sample 1
    # n2                : scalar numeric value defining the number of observations in Sample 2
    # B                 : scalar numeric value indicating the number of replications to use
    # tobit_val         : scalar numeric value indicating the tobit threshold value where censoring is occuring
    # a                 : scalar numeric value indicating the truncation value
    # alpha             : scalar numeric value used to set the Type I probability for the non-inferiority test. Error rate should be alpha/2
    # silent             : logical indicator of whether to print results while running the function. default is FALSE

    requireNamespace("future.apply")
    requireNamespace("tcensReg")
    requireNamespace("truncreg")
    requireNamespace("tidyverse")
    #define the mean for population 2
    mu2 <- mu1 + non_inf_margin

    ls_dt <- future.apply::future_lapply(seq_len(B), function(num_reps){ #looping the function over the number of replicates
        y_star_1  <- tcensReg::rtnorm(n = n1, mu = mu1, sd = sd, a = a)
        y_star_2  <- tcensReg::rtnorm(n = n2, mu = mu2, sd = sd, a = a)
        y_star    <- c(y_star_1, y_star_2)
        y_dl      <- cens_method(y_star, method = "DL", tobit_val)
        y_dl_half <- cens_method(y_star, method = "DL_half", tobit_val)
        y_tobit   <- cens_method(y_star, method = "Tobit", tobit_val)
        cens_ind  <- ifelse(y_tobit == tobit_val, 1, 0)
        group <- c(rep(0, n1), rep(1, n2))
        data.frame(y_star, y_dl, y_dl_half, y_tobit, cens_ind, group)
    }, future.seed=rand_seed)
    dt <- array(unlist(ls_dt), dim = c(n1+n2, 6, B),
                dimnames = list(NULL, c("y_star","y_dl", "y_dl_half", "y_tobit", "cens_ind", "group"), NULL))
    #dt is a large array containing all of the different vector values
    #dimensions are (n1+n2) x dt_vars x num_reps
    #dimension one is the number of observations drawn in the total sample which is equal to n1+n2
    #dimension two is the number of variables in the data.frame which includes the outcome variable y for each method, censoring indicator, and group variable
    #dimension three is the number of different mu combinations. This is equal to length(mu1_vec)
    #dimension four is the number of standard deviation values
    #dimension five is the total number of replicates

    #### Calculate Censoring Percentage #####
    #calculate the percent censoring within each group over the number of mu combinations, sd combinations, and number of replicates
    prop_cens_rep <- apply(dt, 3, function(X){
        df <- data.frame(X)
        cens_tbl <- round(prop.table(table(df$cens_ind, df$group), 2), 4)
        return(cens_tbl)
    })

    #averaging the censoring results over the number of replicates
    prop_cens <- apply(prop_cens_rep, 1, mean)[c(2, 4)] # extracting the cens_ind=1 proportions
    names(prop_cens) <- c("Group_1", "Group_2")
    percent_cens <- data.frame(cbind(prop_cens, mu = c(mu1, mu2),
                                     non_inf = non_inf_margin, sd = sd))
    percent_cens$group <- row.names(percent_cens)
    percent_cens_wide <- pivot_wider(percent_cens, names_from=group, values_from=c(mu, prop_cens))
    percent_cens_wide$n1 <- n1
    percent_cens_wide$n2 <- n2
    if(silent==FALSE){
        #calculating the censoring percentage for each scenario
        writeLines("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        writeLines("Censoring Percentage")
        writeLines("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        print(round(percent_cens, 4))
        writeLines("")
        writeLines("")
    }


    ##### Estimating the Parameters for each of the Six Methods
    method_names <- c("Uncens_NT", "GS", "DL", "DL_half", "Tobit", "tcensReg")
    param_rep <- future.apply::future_apply(dt, 3, function(dt_X){
        df <- data.frame(dt_X)
        #fitting the models
        comp_mod <- lm(y_star ~ group, data = df)
        dl_mod <- lm(y_dl ~ group, data = df)
        dl_half_mod <- lm(y_dl_half ~ group, data = df)
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
    #these rep arrays are 12 x B
    #dim1 represent each of the 6 methods of analysis x 2 parameters (delta and sigma)
    #dim2 represents the number of replications
    diff_rep <- param_rep[1:6, ]
    row.names(diff_rep) <- method_names
    sd_rep <- param_rep[7:12, ]
    row.names(sd_rep) <- method_names

    #####################################
    #####Non-Inferiority Testing#########
    #####################################

    sd_rep_adj <- ((n1+n2)/(n1+n2-2)) * sd_rep #applying an adjustment since the mle was biased
    lb_ci <- diff_rep - qnorm(1 - alpha)*(sd_rep_adj*sqrt(1/n1 + 1/n2)) #calculating the lower_bound confidence interval
    reject_avg <- apply(lb_ci>non_inf_margin, 1, sum)/B #if the lb confidence interval is above the non_inf margin then we reject
    reject_df <-data.frame(obs_reject=reject_avg, true_alpha=alpha)
    reject_results <- data.frame(percent_cens_wide, reject_df)
    if(silent==FALSE){
        writeLines("=======================================")
        writeLines("Rejecting the Null Hyptothesis of Non-Inferiority:")
        writeLines("=======================================")
        print(round(reject_results, 4))
    }
    return(reject_results)
}

#testing to make sure this works
cens_diff_sim_noninf(rand_seed = 032020, mu1 = 0.9,
                     non_inf_margin = -0.15, sd = 0.45,
                     n1 = 10, n2 = 10, B = 100, tobit_val = 0.61, a = 0,
                     alpha = 0.05, silent=TRUE)

tictoc::tic()
tnorm_results_noninf_list <- pbapply::pblapply(samp_sizes, function(x){
    result <- cens_diff_sim_noninf(rand_seed = 032020, mu1 = 0.9,
                         non_inf_margin = -0.15, sd = 0.45,
                         n1 = x, n2 = x, B = 10000, tobit_val = 0.61, a = 0,
                         alpha = 0.05, silent=TRUE)
})
tictoc::toc()

names(tnorm_results_noninf_list) <- paste0("n_", samp_sizes)
tnorm_results_noninf <- do.call(rbind, tnorm_results_noninf_list)
