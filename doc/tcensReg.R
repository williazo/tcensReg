## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("williazo/tcensReg")

## ------------------------------------------------------------------------
library(msm) #we will use this package to generate random values from the truncated normal distribution
mu <- 0.5
sigma <- 0.5
a <- 0

y_star <- msm::rtnorm(n = 1000, mean = mu, sd = sigma, lower = a)
range(y_star) #note that the lowerbound will always be non-negative

## ------------------------------------------------------------------------
nu <- 0.25
y <- ifelse(y_star <= nu, nu, y_star)
sum(y == nu)/length(y) #calculating the number of censored observations
dt <- data.frame(y_star, y) #collecting the uncensored and censored data together


## ----echo=FALSE, warning=FALSE, message = FALSE, out.width="40%", fig.align="center"----
library(ggplot2)
library(viridis)
ggplot(data = dt, aes(x = y_star, y = ..density..))+
  geom_histogram(binwidth = 0.25, col = viridis(1),fill = viridis(1), alpha = 0.6)+
  stat_density(bw = 0.5, col = viridis(1), fill = viridis(1), alpha = 0.3)+
  scale_x_continuous(breaks = seq(0, 2, 0.25))+
  ylab("Density")+
  xlab(expression(Y^"*"))+
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid.major = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(size = 10), legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),  legend.key.width = unit(15, units = "mm"), strip.text = element_text(size = 20),
        axis.title = element_text(size = 20))

## ----echo=FALSE, warning=FALSE, message = FALSE, out.width="40%", fig.align="center"----
ggplot(data = dt, aes(x = y, y = ..density..))+
  geom_histogram(binwidth = 0.25, col = viridis(1, begin = 0.5),fill = viridis(1, begin = 0.5), alpha = 0.6)+
  stat_density(bw = 0.25, col = viridis(1, begin = 0.5), fill = viridis(1, begin = 0.5), alpha = 0.3)+
  scale_x_continuous(breaks = seq(0, 2, 0.25))+
  ylab("Density")+
  xlab(expression(Y))+
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.grid.major = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(size = 10), legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),  legend.key.width = unit(15, units = "mm"), strip.text = element_text(size = 20),
        axis.title = element_text(size = 20))

## ----message = FALSE-----------------------------------------------------
library(tcensReg)  #loading the package into the current environment
tcensReg(y ~ 1, data = dt, a = 0, v = 0.25)

## ----message = FALSE-----------------------------------------------------
output <- tcensReg(y ~ 1, data = dt, a = a, v = nu)
lm_output <- lm(y ~ 1, data = dt) #running OLS model for comparison
cens_output <- tcensReg(y ~ 1, data = dt, v = nu) #censored only model, i.e., Tobit model

tcensReg_est <- output$theta #extracting the point estimates
tcensReg_est[2] <- exp(tcensReg_est[2]) #exponentiating the estimate of log_sigma to estimate sigma

lm_est <- c(coef(lm_output), summary(lm_output)$sigma)

cens_est <- cens_output$theta
cens_est[2] <- exp(cens_est[2])

results_df <- data.frame(rbind(c(mu, sigma), t(tcensReg_est), lm_est, t(cens_est)))
names(results_df) <- c("mu", "sigma")
row.names(results_df) <- c("Truth", "tcensReg", "Normal MLE", "Tobit")
results_df$mu_bias <- abs(results_df$mu - mu)
results_df$sigma_bias <- abs(results_df$sigma - sigma)

knitr::kable(results_df, format = "markdown", digits = 4)

