library(tidyverse)
library(MASS)
library(np)
library(moments)
library(nortest)

########
# This script accompanies the empirical example in the publication:
# "Dealing with Regression Models' Endogeneity by Means of an Adjusted Estimator for the Gaussian Copula Approach".
# Please note that the data used here is simulated; therefore, the results produced by this script will differ from those reported in the article.
# The simulated data is provided to allow the script to run and to enable users to explore the functions and procedures
# discussed in the article.
# For access to the original dataset, please contact the authors of the article:
# Park, S., & Gupta, S. (2012). Handling Endogenous Regressors by Joint Estimation Using Copulas. Marketing Science, 31(4), 567-586.
########

# Functions ---------------------------------------------------------------
# Proposed ECDF
cdf_est_f4 = function(x){
  n <- length(x)
  1/(2*n) + (n-1)/(n)* ecdf(x)(x)
}
# Calculate marginal copulas over each dummy category
get_copula <- function(endog, categories = NULL){
  if (is.null(categories)) {
    copula <- qnorm(cdf_est_f4(endog))  
  }else{
    copula <- matrix(0, nrow = length(categories[[1]]))
    for (c in 1:length(categories)) {
      copula[categories[[c]] == 1] <- qnorm(cdf_est_f4(endog[categories[[c]] == 1]))
    }
  }
  return(copula)
}
# Reduce copula dimension if exogenous continuous variable is in the model
cop_dim_reduce <- function(endog_cop, exog_cop = NULL, Z = NULL){
  # Creating copula separate for a value of Z
  if (!is.null(Z)) {
    # Select only observations in category of Z
    endog_cop_sub <- lapply(endog_cop, FUN = function(x) x[Z==1,])
    # Only create W copula if a W exists
    if (!is.null(exog_cop) ) {
      exog_cop_sub <- lapply(exog_cop, FUN = function(x) x[Z==1,] )  
    }else{exog_cop_sub <- list(NULL)}
    
  }else{ # Create copula for all observations
    endog_cop_sub <- endog_cop
    # Only create W copula if a W exists
    if (!is.null(exog_cop) ) {
      exog_cop_sub <- exog_cop
    }else{exog_cop_sub <- list(NULL)}
  }
  
  # Matrix of marginal copulas
  C_P_W_vec <- cbind(do.call(cbind, endog_cop_sub), do.call(cbind, exog_cop_sub))
  I_dp <- diag(length(endog_cop_sub))
  # Only dw_dp_0 if we have exogenous variables with a copula structure
  if (!is.null(exog_cop)) {
    dw_dp_0 <- matrix(0, nrow = length(exog_cop_sub), ncol = length(endog_cop_sub))  
  }else{ 
    dw_dp_0 <- NULL
  }
  
  # Prepare vector for copulas
  C_P_W <- C_P_W_vec %*% ginv(var(C_P_W_vec)) %*% rbind(I_dp, dw_dp_0)
  # Naming
  colnames(C_P_W) <- names(endog_cop)
  return(C_P_W)
}
# Function for bootstrap inference
boot_reg <- function(dat, B, ols_formula){
  # Regression with original data
  org_reg <- lm(ols_formula, data = dat)
  # Matrices to hold bootstrap results 
  boot_coef <- boot_tstat_h0 <- matrix(NA, nrow = B, ncol = length(coef(org_reg)),
                                       dimnames = list(NULL, names(coef(org_reg))))
  # Running B bootstrap resamples
  for (b in 1:B) {
    boot_index <- sample(1:nrow(dat), replace = TRUE)
    boot_reg <- lm(ols_formula, data = dat[boot_index,])
    boot_coef[b,] <- coef(boot_reg)
    boot_tstat_h0[b,] <- (coef(boot_reg) - coef(org_reg))/summary(boot_reg)$coefficients[,c("Std. Error")]
  }
  # Boot standard error
  boot_std_error <- apply(boot_coef, MARGIN = 2, FUN = function(x) sqrt(var(x)))
  # t-stat using bootstrapped standard errors
  boot_t_stats <- coef(org_reg)/boot_std_error
  # P-values from evaluating t-stat using bootstrapped standard errors against t-distribution
  boot_p_val_from_boot_std_error <- boot_t_stats*NA
  for (c in 1:length(boot_p_val_from_boot_std_error)) {
    boot_p_val_from_boot_std_error[c] <- round(2 * (1 - pt(abs(boot_t_stats[c]), summary(org_reg)$df[2])), digits = 5)  
  }
  
  # P-value based on bootstrapped t-statistics
  # See e.g. "introduction to the bootstrap" by Efron and Tibshirani 1993 (p. 224), or
  # "Microeconometrics" by Cameron and Trivedi 2005 (p. 363).
  boot_percentile_p_val <- matrix(NA, nrow = 1, ncol = length(coef(org_reg)),
                                  dimnames = list(NULL, names(coef(org_reg))))
  
  for (c in 1:ncol(boot_percentile_p_val)) {
    tmp_org_res <- summary(org_reg)$coefficients[c, c("t value")]
    tmp_boot_res <- boot_tstat_h0[,c]
    # Fraction of bootstrap estimates more extreme than original estimates  
    boot_percentile_p_val[,c] <- (sum(tmp_boot_res > abs(tmp_org_res)) + sum(tmp_boot_res<= -abs(tmp_org_res)))/B 
  }
  # Percentile confidence intervals
  boot_percentile_conf <- matrix(NA, nrow = 1, ncol = length(coef(org_reg)),
                                 dimnames = list(NULL, names(coef(org_reg))))
  for (c in 1:ncol(boot_percentile_conf)) {
    tmp_org_res <- summary(org_reg)$coefficients[c, c("t value")]
    tmp_boot_res <- boot_tstat_h0[,c]
    # Fraction of bootstrap estimates more extreme than original estimates  
    boot_percentile_conf[,c] <-  paste0("(",paste(round(quantile(boot_coef[,c], prob = c(0.025, 0.975)),digits = 3), collapse = ";"),")")
  }
  
  return(list(boot_std_error = boot_std_error,
              boot_t_stats = boot_t_stats,
              boot_p_val_from_boot_std_error = boot_p_val_from_boot_std_error,
              boot_percentile_p_val = boot_percentile_p_val,
              boot_percentile_conf = boot_percentile_conf))
  
}
boot_reg_PG_scale <- function(dat, B, ols_formula){
  bws_select <- npudistbw(dat$Price, bwmethod = "cv.cdf", ckertype = "epanechnikov", ckerorder=2)
  bws_select$bw <- 0.0225
  est_cdf <- npudist(bws = bws_select, edat = dat$Price, ckertype = "epanechnikov")
  dat$copula1 <- qnorm(fitted(est_cdf))
  # Regression with original data
  pg_no_cop <- lm(YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4, data = dat)
  varres <- sqrt(var(pg_no_cop$residuals))
  org_reg <- lm(YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4 + I(copula1*varres), data = dat) 
  # Matrices to hold bootstrap results 
  boot_coef <- boot_tstat_h0 <- matrix(NA, nrow = B, ncol = length(coef(org_reg)),
                                       dimnames = list(NULL, names(coef(org_reg))))
  # Running B bootstrap resamples
  for (b in 1:B) {
    boot_dat <- dat[sample(1:nrow(dat), replace = TRUE),]
    boot_no_cop <- lm(YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4, data = boot_dat)
    boot_varres <- sqrt(var(boot_no_cop$residuals))
    
    est_cdf <- npudist(bws = bws_select, edat = boot_dat$Price, ckertype = "epanechnikov")
    boot_dat$copula1 <- qnorm(fitted(est_cdf))
    
    boot_reg <- lm(YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4 + I(copula1*boot_varres), data = boot_dat) 
    boot_coef[b,] <- coef(boot_reg)
    boot_tstat_h0[b,] <- (coef(boot_reg) - coef(org_reg))/summary(boot_reg)$coefficients[,c("Std. Error")]
  }
  # Boot standard error
  boot_std_error <- apply(boot_coef, MARGIN = 2, FUN = function(x) sqrt(var(x)))
  # t-stat using bootstrapped standard errors
  boot_t_stats <- coef(org_reg)/boot_std_error
  # P-values from evaluating t-stat using bootstrapped standard errors against t-distribution
  boot_p_val_from_boot_std_error <- boot_t_stats*NA
  for (c in 1:length(boot_p_val_from_boot_std_error)) {
    boot_p_val_from_boot_std_error[c] <- round(2 * (1 - pt(abs(boot_t_stats[c]), summary(org_reg)$df[2])), digits = 5)  
  }
  
  # P-value based on bootstrapped t-statistics
  # See e.g. "introduction to the bootstrap" by Efron and Tibshirani 1993 (p. 224), or
  # "Microeconometrics" by Cameron and Trivedi 2005 (p. 363).
  boot_percentile_p_val <- matrix(NA, nrow = 1, ncol = length(coef(org_reg)),
                                  dimnames = list(NULL, names(coef(org_reg))))
  
  for (c in 1:ncol(boot_percentile_p_val)) {
    tmp_org_res <- summary(org_reg)$coefficients[c, c("t value")]
    tmp_boot_res <- boot_tstat_h0[,c]
    # Fraction of bootstrap estimates more extreme than original estimates  
    boot_percentile_p_val[,c] <- (sum(tmp_boot_res > abs(tmp_org_res)) + sum(tmp_boot_res<= -abs(tmp_org_res)) )/B 
  }
  # Percentile confidence intervals
  boot_percentile_conf <- matrix(NA, nrow = 1, ncol = length(coef(org_reg)),
                                 dimnames = list(NULL, names(coef(org_reg))))
  for (c in 1:ncol(boot_percentile_conf)) {
    tmp_org_res <- summary(org_reg)$coefficients[c, c("t value")]
    tmp_boot_res <- boot_tstat_h0[,c]
    # Fraction of bootstrap estimates more extreme than original estimates  
    boot_percentile_conf[,c] <-  paste0("(",paste(round(quantile(boot_coef[,c], prob = c(0.025, 0.975)),digits = 3), collapse = ";"),")")
  }
  
  return(list(boot_std_error = boot_std_error,
              boot_t_stats = boot_t_stats,
              boot_p_val_from_boot_std_error = boot_p_val_from_boot_std_error,
              boot_percentile_p_val = boot_percentile_p_val,
              boot_percentile_conf = boot_percentile_conf))
  
}
boot_reg_PG <- function(dat, B, ols_formula){
  bws_select <- npudistbw(dat$Price, bwmethod = "cv.cdf", ckertype = "epanechnikov", ckerorder=2)
  bws_select$bw <- 0.0225
  est_cdf <- npudist(bws = bws_select, edat = dat$Price, ckertype = "epanechnikov")
  dat$copula1 <- qnorm(fitted(est_cdf))
  # Original regression
  org_reg <- lm(YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4 + copula1, data = dat) 
  # Matrices to hold bootstrap results 
  boot_coef <- boot_tstat_h0 <- matrix(NA, nrow = B, ncol = length(coef(org_reg)),
                                       dimnames = list(NULL, names(coef(org_reg))))
  bws_select <- npudistbw(dat$Price, bwmethod = "cv.cdf", ckertype = "epanechnikov", ckerorder=2)
  bws_select$bw <- 0.0225
  # Running B bootstrap resamples
  for (b in 1:B) {
    boot_dat <- dat[sample(1:nrow(dat), replace = TRUE),]
    
    est_cdf <- npudist(bws = bws_select, edat = boot_dat$Price, ckertype = "epanechnikov")
    boot_dat$copula1 <- qnorm(fitted(est_cdf))
    
    boot_reg <- lm(YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4 + copula1, data = boot_dat) 
    boot_coef[b,] <- coef(boot_reg)
    boot_tstat_h0[b,] <- (coef(boot_reg) - coef(org_reg))/summary(boot_reg)$coefficients[,c("Std. Error")]
  }
  # Boot standard error
  boot_std_error <- apply(boot_coef, MARGIN = 2, FUN = function(x) sqrt(var(x)))
  # t-stat using bootstrapped standard errors
  boot_t_stats <- coef(org_reg)/boot_std_error
  # P-values from evaluating t-stat using bootstrapped standard errors against t-distribution
  boot_p_val_from_boot_std_error <- boot_t_stats*NA
  for (c in 1:length(boot_p_val_from_boot_std_error)) {
    boot_p_val_from_boot_std_error[c] <- round(2 * (1 - pt(abs(boot_t_stats[c]), summary(org_reg)$df[2])), digits = 5)  
  }
  
  # P-value based on bootstrapped t-statistics
  # See e.g. "introduction to the bootstrap" by Efron and Tibshirani 1993 (p. 224), or
  # "Microeconometrics" by Cameron and Trivedi 2005 (p. 363).
  boot_percentile_p_val <- matrix(NA, nrow = 1, ncol = length(coef(org_reg)),
                                  dimnames = list(NULL, names(coef(org_reg))))
  
  for (c in 1:ncol(boot_percentile_p_val)) {
    tmp_org_res <- summary(org_reg)$coefficients[c, c("t value")]
    tmp_boot_res <- boot_tstat_h0[,c]
    # Fraction of bootstrap estimates more extreme than original estimates  
    boot_percentile_p_val[,c] <- (sum(tmp_boot_res > abs(tmp_org_res)) + sum(tmp_boot_res<= -abs(tmp_org_res)) )/B 
  }
  # Percentile confidence intervals
  boot_percentile_conf <- matrix(NA, nrow = 1, ncol = length(coef(org_reg)),
                                 dimnames = list(NULL, names(coef(org_reg))))
  for (c in 1:ncol(boot_percentile_conf)) {
    tmp_org_res <- summary(org_reg)$coefficients[c, c("t value")]
    tmp_boot_res <- boot_tstat_h0[,c]
    # Fraction of bootstrap estimates more extreme than original estimates  
    boot_percentile_conf[,c] <-  paste0("(",paste(round(quantile(boot_coef[,c], prob = c(0.025, 0.975)),digits = 3), collapse = ";"),")")
  }
  
  return(list(boot_std_error = boot_std_error,
              boot_t_stats = boot_t_stats,
              boot_p_val_from_boot_std_error = boot_p_val_from_boot_std_error,
              boot_percentile_p_val = boot_percentile_p_val,
              boot_percentile_conf = boot_percentile_conf))
  
}
boot_proposed_price_endog <- function(dat, B, ols_formula, coef_joint_test){
  # Regression with original data
  org_reg <- lm(ols_formula, data = dat)
  org_coef <- org_reg$coefficients[names(org_reg$coefficients) %in% coef_joint_test]
  beta_vcov <- vcov(org_reg)
  beta_vcov <- beta_vcov[rownames(beta_vcov) %in% coef_joint_test , colnames(beta_vcov) %in% coef_joint_test]
  org_wald_test <- t(org_coef) %*% ginv(beta_vcov) %*% (org_coef)
  # Matrices to hold bootstrap results 
  boot_coef <- boot_tstat_h0 <- matrix(NA, nrow = B, ncol = length(coef(org_reg)),
                                       dimnames = list(NULL, names(coef(org_reg))))
  boot_wald_test <- matrix(NA, nrow = B)
  # Running B bootstrap resamples
  for (b in 1:B) {
    boot_dat <- dat[sample(1:nrow(dat), replace = TRUE),]
    # List of dummy categories
    categories <- list(boot_dat$Q1, boot_dat$Q2, boot_dat$Q3, boot_dat$Q4)
    # Create marginal copula terms within each category
    boot_dat$cop_Price <- get_copula(endog = boot_dat$Price, categories)
    boot_dat$cop_Promo <- get_copula(endog = boot_dat$Promo, categories)
    boot_dat$cop_Disp <- get_copula(endog = boot_dat$Disp, categories)
    # Create C_P matrix
    endog_cop <- list(copula1 = boot_dat$cop_Price)
    exo_cop <- list(copula1 = boot_dat$cop_Promo, copula2 = boot_dat$cop_Disp)
    C_P_W <- matrix(NA, nrow = nrow(boot_dat), ncol = length(endog_cop))
    colnames(C_P_W) <- names(endog_cop)
    C_P_W[boot_dat$Q1==1,] <- cop_dim_reduce(endog_cop, exo_cop, Z = boot_dat$Q1)
    C_P_W[boot_dat$Q2==1,] <- cop_dim_reduce(endog_cop, exo_cop, Z = boot_dat$Q2)
    C_P_W[boot_dat$Q3==1,] <- cop_dim_reduce(endog_cop, exo_cop, Z = boot_dat$Q3)
    C_P_W[boot_dat$Q4==1,] <- cop_dim_reduce(endog_cop, exo_cop, Z = boot_dat$Q4)
    boot_dat$C_P_W <- C_P_W
    
    boot_reg <- lm(ols_formula, data = boot_dat)
    boot_coef[b,] <- coef(boot_reg)
    boot_tstat_h0[b,] <- (coef(boot_reg) - coef(org_reg))/summary(boot_reg)$coefficients[,c("Std. Error")]
    # Test that all copula terms are zero
    subset_boot_coef <- boot_reg$coefficients[names(boot_reg$coefficients) %in% coef_joint_test]
    boot_beta_vcov <- vcov(boot_reg)
    boot_beta_vcov <- boot_beta_vcov[rownames(boot_beta_vcov) %in% coef_joint_test , colnames(boot_beta_vcov) %in% coef_joint_test]
    boot_wald_test[b,] <- t(subset_boot_coef - org_coef) %*% ginv(boot_beta_vcov) %*% (subset_boot_coef - org_coef)
    
  }
  # Boot standard error
  boot_std_error <- apply(boot_coef, MARGIN = 2, FUN = function(x) sqrt(var(x)))
  # t-stat using bootstrapped standard errors
  boot_t_stats <- coef(org_reg)/boot_std_error
  # P-values from evaluating t-stat using bootstrapped standard errors against t-distribution
  boot_p_val_from_boot_std_error <- boot_t_stats*NA
  for (c in 1:length(boot_p_val_from_boot_std_error)) {
    boot_p_val_from_boot_std_error[c] <- round(2 * (1 - pt(abs(boot_t_stats[c]), summary(org_reg)$df[2])), digits = 5)  
  }
  
  # P-value based on bootstrapped t-statistics
  # See e.g. "introduction to the bootstrap" by Efron and Tibshirani 1993 (p. 224), or
  # "Microeconometrics" by Cameron and Trivedi 2005 (p. 363).
  boot_percentile_p_val <- matrix(NA, nrow = 1, ncol = length(coef(org_reg)),
                                  dimnames = list(NULL, names(coef(org_reg))))
  
  for (c in 1:ncol(boot_percentile_p_val)) {
    tmp_org_res <- summary(org_reg)$coefficients[c, c("t value")]
    tmp_boot_res <- boot_tstat_h0[,c]
    # Fraction of bootstrap estimates more extreme than original estimates  
    boot_percentile_p_val[,c] <- (sum(tmp_boot_res > abs(tmp_org_res)) + sum(tmp_boot_res<= -abs(tmp_org_res)) )/B 
  }
  # Bootstrap confidence intervals
  boot_percentile_conf <- matrix(NA, nrow = 1, ncol = length(coef(org_reg)),
                                 dimnames = list(NULL, names(coef(org_reg))))
  for (c in 1:ncol(boot_percentile_conf)) {
    tmp_org_res <- summary(org_reg)$coefficients[c, c("t value")]
    tmp_boot_res <- boot_tstat_h0[,c]
    # Fraction of bootstrap estimates more extreme than original estimates  
    boot_percentile_conf[,c] <-  paste0("(",paste(round(quantile(boot_coef[,c], prob = c(0.025, 0.975)),digits = 3), collapse = ";"),")")
  }
  
  # P-value for bootstrapped Wald test
  wald_test_p_val <- (sum(org_wald_test[1,1] < boot_wald_test[,1]))/B
  
  
  
  return(list(boot_std_error = boot_std_error,
              boot_t_stats = boot_t_stats,
              boot_p_val_from_boot_std_error = boot_p_val_from_boot_std_error,
              boot_percentile_p_val = boot_percentile_p_val,
              boot_percentile_conf = boot_percentile_conf,
              wald_test_p_val = wald_test_p_val))
  
}
boot_proposed_all_endog <- function(dat, B, ols_formula, coef_joint_test){
  # Regression with original data
  org_reg <- lm(ols_formula, data = dat)
  org_coef <- org_reg$coefficients[names(org_reg$coefficients) %in% coef_joint_test]
  beta_vcov <- vcov(org_reg)
  #beta_vcov <- vcovHC(org_reg)
  beta_vcov <- beta_vcov[rownames(beta_vcov) %in% coef_joint_test , colnames(beta_vcov) %in% coef_joint_test]
  org_wald_test <- t(org_coef) %*% ginv(beta_vcov) %*% (org_coef)
  # Matrices to hold bootstrap results 
  boot_coef <- boot_tstat_h0 <- matrix(NA, nrow = B, ncol = length(coef(org_reg)),
                                       dimnames = list(NULL, names(coef(org_reg))))
  boot_wald_test <- matrix(NA, nrow = B)
  # Running B bootstrap resamples
  for (b in 1:B) {
    boot_dat <- dat[sample(1:nrow(dat), replace = TRUE),]
    # List of dummy categories
    categories <- list(boot_dat$Q1, boot_dat$Q2, boot_dat$Q3, boot_dat$Q4)
    # Create marginal copula terms within each category
    boot_dat$cop_Price <- get_copula(endog = boot_dat$Price, categories)
    boot_dat$cop_Promo <- get_copula(endog = boot_dat$Promo, categories)
    boot_dat$cop_Disp <- get_copula(endog = boot_dat$Disp, categories)
    # Create C_P matrix
    endog_cop <- list(copula1 = boot_dat$cop_Price, copula2 = boot_dat$cop_Promo, copula3 = boot_dat$cop_Disp)
    C_P_W <- matrix(NA, nrow = nrow(boot_dat), ncol = length(endog_cop))
    colnames(C_P_W) <- names(endog_cop)
    C_P_W[boot_dat$Q1==1,] <- cop_dim_reduce(endog_cop, Z = boot_dat$Q1)
    C_P_W[boot_dat$Q2==1,] <- cop_dim_reduce(endog_cop, Z = boot_dat$Q2)
    C_P_W[boot_dat$Q3==1,] <- cop_dim_reduce(endog_cop, Z = boot_dat$Q3)
    C_P_W[boot_dat$Q4==1,] <- cop_dim_reduce(endog_cop, Z = boot_dat$Q4)
    # Bind to boot_dat
    boot_dat <- cbind(boot_dat, C_P_W)
    # Boot regression
    boot_reg <- lm(ols_formula, data = boot_dat)
    boot_coef[b,] <- coef(boot_reg)
    boot_tstat_h0[b,] <- (coef(boot_reg) - coef(org_reg))/summary(boot_reg)$coefficients[,c("Std. Error")]
    # Test that all copula terms are zero
    subset_boot_coef <- boot_reg$coefficients[names(boot_reg$coefficients) %in% coef_joint_test]
    boot_beta_vcov <- vcov(boot_reg)
    #boot_beta_vcov <- vcovHC(boot_reg)
    
    boot_beta_vcov <- boot_beta_vcov[rownames(boot_beta_vcov) %in% coef_joint_test , colnames(boot_beta_vcov) %in% coef_joint_test]
    boot_wald_test[b,] <- t(subset_boot_coef - org_coef) %*% ginv(boot_beta_vcov) %*% (subset_boot_coef - org_coef)
    
  }
  # Boot standard error
  boot_std_error <- apply(boot_coef, MARGIN = 2, FUN = function(x) sqrt(var(x)))
  # t-stat using bootstrapped standard errors
  boot_t_stats <- coef(org_reg)/boot_std_error
  # P-values from evaluating t-stat using bootstrapped standard errors against t-distribution
  boot_p_val_from_boot_std_error <- boot_t_stats*NA
  for (c in 1:length(boot_p_val_from_boot_std_error)) {
    boot_p_val_from_boot_std_error[c] <- round(2 * (1 - pt(abs(boot_t_stats[c]), summary(org_reg)$df[2])), digits = 5)  
  }
  
  # P-value based on bootstrapped t-statistics
  # See e.g. "introduction to the bootstrap" by Efron and Tibshirani 1993 (p. 224), or
  # "Microeconometrics" by Cameron and Trivedi 2005 (p. 363).
  boot_percentile_p_val <- matrix(NA, nrow = 1, ncol = length(coef(org_reg)),
                                  dimnames = list(NULL, names(coef(org_reg))))
  
  for (c in 1:ncol(boot_percentile_p_val)) {
    tmp_org_res <- summary(org_reg)$coefficients[c, c("t value")]
    tmp_boot_res <- boot_tstat_h0[,c]
    # Fraction of bootstrap estimates more extreme than original estimates  
    boot_percentile_p_val[,c] <- (sum(tmp_boot_res > abs(tmp_org_res)) + sum(tmp_boot_res<= -abs(tmp_org_res)) )/B 
  }
  # Bootstrap confidence intervals
  boot_percentile_conf <- matrix(NA, nrow = 1, ncol = length(coef(org_reg)),
                                 dimnames = list(NULL, names(coef(org_reg))))
  for (c in 1:ncol(boot_percentile_conf)) {
    tmp_org_res <- summary(org_reg)$coefficients[c, c("t value")]
    tmp_boot_res <- boot_tstat_h0[,c]
    # Fraction of bootstrap estimates more extreme than original estimates  
    boot_percentile_conf[,c] <-  paste0("(",paste(round(quantile(boot_coef[,c], prob = c(0.025, 0.975)),digits = 3), collapse = ";"),")")
  }
  # P-value for bootstrapped Wald test
  wald_test_p_val <- (sum(org_wald_test[1,1] < boot_wald_test[,1]))/B
  
  return(list(boot_std_error = boot_std_error,
              boot_t_stats = boot_t_stats,
              boot_p_val_from_boot_std_error = boot_p_val_from_boot_std_error,
              boot_percentile_p_val = boot_percentile_p_val,
              boot_percentile_conf = boot_percentile_conf,
              wald_test_p_val = wald_test_p_val))
  
}
# Checking non-normality of log(Price), Promo and Display -----------------

####
# Store 1
####
# Read data
pg_dat <- readxl::read_excel(path = "simulated_data.xlsx", sheet = 1)
# Skewness 
moments::skewness(pg_dat$Price)
moments::skewness(pg_dat$Promo)
moments::skewness(pg_dat$Disp)
# Plot densities
ggplot(pg_dat, aes(x=Price)) +
  geom_density(alpha=0.4, fill="red") +
  stat_function(fun = dnorm,
                args = list(mean = mean(pg_dat$Price),
                            sd = sd(pg_dat$Price)), 
                color="blue") +
  xlim(-1,1) +
  labs(y = "Density of Price", x = "Price") + 
  theme_classic()


ggplot(pg_dat, aes(x=Promo)) +
  geom_density(alpha=0.4, fill="red") +
  stat_function(fun = dnorm,
                args = list(mean = mean(pg_dat$Promo),
                            sd = sd(pg_dat$Promo)), 
                color="blue") +
  xlim(-2.5,2.5) +
  labs(y = "Density of Promotion", x = "Promotion") + 
  theme_classic()


ggplot(pg_dat, aes(x=Disp)) +
  geom_density(alpha=0.4, fill="red") +
  stat_function(fun = dnorm,
                args = list(mean = mean(pg_dat$Disp),
                            sd = sd(pg_dat$Disp)), 
                color="blue") +
  xlim(-1.5,1.5) +
  labs(y = "Density of Display", x = "Display") +
  theme_classic()


# Anderson-Darling test for normality
nortest::ad.test(pg_dat$Price)$p.value
nortest::ad.test(pg_dat$Promo)$p.value
nortest::ad.test(pg_dat$Disp)$p.value
# Cramer-von Mises test for normality
nortest::cvm.test(pg_dat$Price)$p.value
nortest::cvm.test(pg_dat$Promo)$p.value
nortest::cvm.test(pg_dat$Disp)$p.value

####
# Store 2
####
# Read data
pg_dat <- readxl::read_excel(path = "simulated_data.xlsx", sheet = 2)
# Skewness 
moments::skewness(pg_dat$Price)
moments::skewness(pg_dat$Promo)
moments::skewness(pg_dat$Disp)
# Plotting
ggplot(pg_dat, aes(x=Price)) +
  geom_density(alpha=0.4, fill="red") +
  stat_function(fun = dnorm,
                args = list(mean = mean(pg_dat$Price),
                            sd = sd(pg_dat$Price)), 
                color="blue") +
  xlim(-1,1) +
  theme_classic()

ggplot(pg_dat, aes(x=Promo)) +
  geom_density(alpha=0.4, fill="red") +
  stat_function(fun = dnorm,
                args = list(mean = mean(pg_dat$Promo),
                            sd = sd(pg_dat$Promo)), 
                color="blue") +
  xlim(-1.5,2.5) +
  theme_classic()

ggplot(pg_dat, aes(x=Disp)) +
  geom_density(alpha=0.4, fill="red") +
  stat_function(fun = dnorm,
                args = list(mean = mean(pg_dat$Disp),
                            sd = sd(pg_dat$Disp)), 
                color="blue") +
  xlim(-0.5,1.5) +
  theme_classic()

# Anderson-Darling test for normality
nortest::ad.test(pg_dat$Price)$p.value
nortest::ad.test(pg_dat$Promo)$p.value
nortest::ad.test(pg_dat$Disp)$p.value
# Cramer-von Mises test for normality
nortest::cvm.test(pg_dat$Price)$p.value
nortest::cvm.test(pg_dat$Promo)$p.value
nortest::cvm.test(pg_dat$Disp)$p.value

# OLS ---------------------------------------------------------------------
####
# Store 1
####
pg_dat <- readxl::read_excel(path = "simulated_data.xlsx", sheet = 1)
# Add quarter 1 dummy
pg_dat[,"Q1"] <- 0
pg_dat[(pg_dat$Q2 == 0 & pg_dat$Q3 == 0) & (pg_dat$Q4 == 0),"Q1"] <- 1
# OLS formula
ols_formula <- "YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4"
# OLS 
reg_ols <- lm(ols_formula, data = pg_dat)
coef_ols <- coef(reg_ols)
print("OLS store 1: Estimated coefficients")
round(coef_ols, digits = 3)
# OLS bootstrap inference
set.seed(1)
boot_res_ols_store1 <- boot_reg(dat = pg_dat, B = 5000, ols_formula = ols_formula)
print("OLS store 1: Bootstrapped standarderrors")
round(boot_res_ols_store1$boot_std_error, digits = 3)
print("OLS store 1: P-values using bootstrapped standarderrors")
round(boot_res_ols_store1$boot_p_val_from_boot_std_error, digits = 3)
print("OLS store 1: Bootstrap percentile confidence intervals")
boot_res_ols_store1$boot_percentile_conf # Percentile confidence intervals

####
# Store 2
####
pg_dat <- readxl::read_excel(path = "simulated_data.xlsx", sheet = 2)
# Add quarter 1 dummy
pg_dat[,"Q1"] <- 0
pg_dat[(pg_dat$Q2 == 0 & pg_dat$Q3 == 0) & (pg_dat$Q4 == 0),"Q1"] <- 1
# OLS formula
ols_formula <- "YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4"
# OLS
reg_ols <- lm(ols_formula, data = pg_dat)
coef_ols <- coef(reg_ols)
print("OLS store 2: Estimated coefficients")
round(coef_ols, digits = 3)
# OLS bootstrap inference
set.seed(1)
boot_res_ols_store2 <- boot_reg(dat = pg_dat, B = 5000, ols_formula = ols_formula)

print("OLS store 2: Bootstrapped standarderrors")
round(boot_res_ols_store2$boot_std_error, digits = 3)
print("OLS store 2: P-values using bootstrapped standarderrors")
round(boot_res_ols_store2$boot_p_val_from_boot_std_error, digits = 3)
print("OLS store 2: Bootstrap percentile confidence intervals")
boot_res_ols_store2$boot_percentile_conf # Percentile confidence intervals

# P&G replication ---------------------------------------------------------
####
# Store 1
####
pg_dat <- readxl::read_excel(path = "simulated_data.xlsx", sheet = 1)
# Add quarter 1 dummy
pg_dat[,"Q1"] <- 0
pg_dat[(pg_dat$Q2 == 0 & pg_dat$Q3 == 0) & (pg_dat$Q4 == 0),"Q1"] <- 1
# With epanechnikov kernel
bws_select <- npudistbw(pg_dat$Price, bwmethod = "cv.cdf", ckertype = "epanechnikov", ckerorder=2)
bws_select$bw <- 0.0225
est_cdf <- npudist(bws = bws_select, edat = pg_dat$Price, ckertype = "epanechnikov")
pg_dat$copula1 <- qnorm(fitted(est_cdf))

# P&G with non-rescaled copula term (Table 2)
pg_store1 <- lm(YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4 + copula1, data = pg_dat) 
print("P&G replication store 1: Estimated coefficients")
round(pg_store1$coefficients, digits = 3)
ols_formula <- "YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4 + copula1"
set.seed(1)
boot_res_pg_store1 <- boot_reg_PG(dat = pg_dat, B = 5000, ols_formula = ols_formula)

print("P&G replication store 1: Bootstrapped standarderrors")
round(boot_res_pg_store1$boot_std_error, digits = 3)
print("P&G replication store 1: P-values using bootstrapped standarderrors")
round(boot_res_pg_store1$boot_p_val_from_boot_std_error, digits = 3)
print("P&G replication store 1: Bootstrap percentile confidence intervals")
boot_res_pg_store1$boot_percentile_conf

# P&G with rescaled copula term (footnote related to P&G replication)
pg_store1 <- lm(YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4, data = pg_dat)
varres <- sqrt(var(pg_store1$residuals))
pg_store1 <- lm(YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4 + I(copula1*varres), data = pg_dat) 
round(pg_store1$coefficients, digits = 3)

ols_formula <- "YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4 + I(copula1*varres)"
set.seed(1)
boot_res_pg_scale_store1 <- boot_reg_PG_scale(dat = pg_dat, B = 5000, ols_formula = ols_formula)
round(boot_res_pg_scale_store1$boot_p_val_from_boot_std_error, digits = 3)

####
# Store 2
####
pg_dat <- readxl::read_excel(path = "simulated_data.xlsx", sheet = 2)
# Add quarter 1 dummy
pg_dat[,"Q1"] <- 0
pg_dat[(pg_dat$Q2 == 0 & pg_dat$Q3 == 0) & (pg_dat$Q4 == 0),"Q1"] <- 1
# With epanechnikov kernel
bws_select <- npudistbw(pg_dat$Price, bwmethod = "cv.cdf", ckertype = "epanechnikov", ckerorder=2)
bws_select$bw <- 0.0225
est_cdf <- npudist(bws = bws_select, edat = pg_dat$Price, ckertype = "epanechnikov")
pg_dat$copula1 <- qnorm(fitted(est_cdf))

# P&G with non-rescaled copula term  (Table 2)
pg_store2 <- lm(YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4 + copula1, data = pg_dat) 
print("P&G replication store 2: Estimated coefficients")
round(pg_store2$coefficients, digits = 3)
ols_formula <- "YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4 + copula1"
set.seed(1)
boot_res_pg_store2 <- boot_reg_PG(dat = pg_dat, B = 5000, ols_formula = ols_formula)

print("P&G replication store 2: Bootstrapped standarderrors")
round(boot_res_pg_store2$boot_std_error, digits = 3)
print("P&G replication store 2: P-values using bootstrapped standarderrors")
round(boot_res_pg_store2$boot_p_val_from_boot_std_error, digits = 3)
print("P&G replication store 2: Bootstrap percentile confidence intervals")
boot_res_pg_store2$boot_percentile_conf

# P&G with rescaled copula term (footnote related to P&G replication)
pg_store2 <- lm(YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4, data = pg_dat)
varres <- sqrt(var(pg_store2$residuals))
pg_store2 <- lm(YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4 + I(copula1*varres), data = pg_dat) 
round(pg_store2$coefficients, digits = 3)

ols_formula <- "YMAT ~ Price + Promo + Disp + Q2 + Q3 + Q4 + I(copula1*varres)"
set.seed(1)
boot_res_pg_scale_store2 <- boot_reg_PG_scale(dat = pg_dat, B = 5000, ols_formula = ols_formula)
round(boot_res_pg_scale_store2$boot_p_val_from_boot_std_error, digits = 3)


# Main paper: price as endogenous and copula structure differ across quarters -------------
####
# Store 1
####
pg_dat <- readxl::read_excel(path = "simulated_data.xlsx", sheet = 1)
# Add quarter 1 dummy
pg_dat[,"Q1"] <- 0
pg_dat[(pg_dat$Q2 == 0 & pg_dat$Q3 == 0) & (pg_dat$Q4 == 0),"Q1"] <- 1
# List of dummy categories
categories <- list(pg_dat$Q1, pg_dat$Q2, pg_dat$Q3, pg_dat$Q4)
# Create marginal copula terms within each category
pg_dat$cop_Price <- get_copula(endog = pg_dat$Price, categories)
pg_dat$cop_Promo <- get_copula(endog = pg_dat$Promo, categories)
pg_dat$cop_Disp <- get_copula(endog = pg_dat$Disp, categories)
# Create C_P matrix
endog_cop <- list(copula1 = pg_dat$cop_Price)
exo_cop <- list(copula1 = pg_dat$cop_Promo, copula2 = pg_dat$cop_Disp)
C_P_W <- matrix(NA, nrow = nrow(pg_dat), ncol = length(endog_cop))
colnames(C_P_W) <- names(endog_cop)
C_P_W[pg_dat$Q1==1,] <- cop_dim_reduce(endog_cop, exo_cop, Z = pg_dat$Q1)
C_P_W[pg_dat$Q2==1,] <- cop_dim_reduce(endog_cop, exo_cop, Z = pg_dat$Q2)
C_P_W[pg_dat$Q3==1,] <- cop_dim_reduce(endog_cop, exo_cop, Z = pg_dat$Q3)
C_P_W[pg_dat$Q4==1,] <- cop_dim_reduce(endog_cop, exo_cop, Z = pg_dat$Q4)
pg_dat$C_P_W <- C_P_W
# Regression with copula
ols_formula <- "YMAT ~  Price + Promo + Disp + Q2 + Q3 + Q4 + (Q1 + Q2 + Q3 + Q4):C_P_W"
price_endog_reg_store1 <- lm(ols_formula, data = pg_dat)
print("Price endog store 1: Estimated coefficients")
round(coef(price_endog_reg_store1),digits = 3)
# Bootstrap results 
coef_joint_test <- c("Q1:C_P_W", "Q2:C_P_W", "Q3:C_P_W", "Q4:C_P_W")
set.seed(1)
price_endog_boot_store1 <- boot_proposed_price_endog(dat = pg_dat, 
                                                     B = 5000, 
                                                     ols_formula = ols_formula, 
                                                     coef_joint_test = coef_joint_test)
print("Price endog store 1: Bootstrapped standarderrors")
round(price_endog_boot_store1$boot_std_error, digits = 3)
print("Price endog store 1: P-values using bootstrapped standarderrors")
round(price_endog_boot_store1$boot_p_val_from_boot_std_error, digits = 4)
print("Price endog store 1: Bootstrap percentile confidence intervals")
price_endog_boot_store1$boot_percentile_conf
print("Price endog store 1: p-value from bootstrapped wald stat")
price_endog_boot_store1$wald_test_p_val


####
# Store 2
####
pg_dat <- readxl::read_excel(path = "simulated_data.xlsx", sheet = 2)
# Add quarter 1 dummy
pg_dat[,"Q1"] <- 0
pg_dat[(pg_dat$Q2 == 0 & pg_dat$Q3 == 0) & (pg_dat$Q4 == 0),"Q1"] <- 1
# List of dummy categories
categories <- list(pg_dat$Q1, pg_dat$Q2, pg_dat$Q3, pg_dat$Q4)
# Create marginal copula terms within each category
pg_dat$cop_Price <- get_copula(endog = pg_dat$Price, categories)
pg_dat$cop_Promo <- get_copula(endog = pg_dat$Promo, categories)
pg_dat$cop_Disp <- get_copula(endog = pg_dat$Disp, categories)
# Create C_P matrix
endog_cop <- list(copula1 = pg_dat$cop_Price)
exo_cop <- list(copula1 = pg_dat$cop_Promo, copula2 = pg_dat$cop_Disp)
C_P_W <- matrix(NA, nrow = nrow(pg_dat), ncol = length(endog_cop))
colnames(C_P_W) <- names(endog_cop)
C_P_W[pg_dat$Q1==1,] <- cop_dim_reduce(endog_cop, exo_cop, Z = pg_dat$Q1)
C_P_W[pg_dat$Q2==1,] <- cop_dim_reduce(endog_cop, exo_cop, Z = pg_dat$Q2)
C_P_W[pg_dat$Q3==1,] <- cop_dim_reduce(endog_cop, exo_cop, Z = pg_dat$Q3)
C_P_W[pg_dat$Q4==1,] <- cop_dim_reduce(endog_cop, exo_cop, Z = pg_dat$Q4)
pg_dat$C_P_W <- C_P_W
# Regression with copula
ols_formula <- "YMAT ~  Price + Promo + Disp + Q2 + Q3 + Q4 + (Q1 + Q2 + Q3 + Q4):C_P_W"
price_endog_reg_store2 <- lm(ols_formula, data = pg_dat)
print("Price endog store 2: Estimated coefficients")
round(coef(price_endog_reg_store2), digits = 3)
# Bootstrap results
coef_joint_test <- c("Q1:C_P_W", "Q2:C_P_W", "Q3:C_P_W", "Q4:C_P_W")
set.seed(1)
price_endog_boot_store2 <- boot_proposed_price_endog(dat = pg_dat, 
                                                     B = 5000, 
                                                     ols_formula = ols_formula,
                                                     coef_joint_test = coef_joint_test)
print("Price endog store 2: Bootstrapped standarderrors")
round(price_endog_boot_store2$boot_std_error, digits = 3)
print("Price endog store 2: P-values using bootstrapped standarderrors")
round(price_endog_boot_store2$boot_p_val_from_boot_std_error, digits = 4)
print("Price endog store 2: Bootstrap percentile confidence intervals")
price_endog_boot_store2$boot_percentile_conf
print("Price endog store 2: p-value from bootstrapped wald stat")
price_endog_boot_store2$wald_test_p_val

# Main paper - price, promotion, display are endogenous and copula structure differ across quarters -------------
####
# Store 1
####
pg_dat <- readxl::read_excel(path = "simulated_data.xlsx", sheet = 1)
# Add quarter 1 dummy
pg_dat[,"Q1"] <- 0
pg_dat[(pg_dat$Q2 == 0 & pg_dat$Q3 == 0) & (pg_dat$Q4 == 0),"Q1"] <- 1
# List of dummy categories
categories <- list(pg_dat$Q1, pg_dat$Q2, pg_dat$Q3, pg_dat$Q4)
# Create marginal copula terms within each category
pg_dat$cop_Price <- get_copula(endog = pg_dat$Price, categories)
pg_dat$cop_Promo <- get_copula(endog = pg_dat$Promo, categories)
pg_dat$cop_Disp <- get_copula(endog = pg_dat$Disp, categories)
# Create C_P matrix
endog_cop <- list(copula1 = pg_dat$cop_Price, copula2 = pg_dat$cop_Promo, copula3 = pg_dat$cop_Disp)
C_P_W <- matrix(NA, nrow = nrow(pg_dat), ncol = length(endog_cop))
colnames(C_P_W) <- names(endog_cop)
C_P_W[pg_dat$Q1==1,] <- cop_dim_reduce(endog_cop, Z = pg_dat$Q1)
C_P_W[pg_dat$Q2==1,] <- cop_dim_reduce(endog_cop, Z = pg_dat$Q2)
C_P_W[pg_dat$Q3==1,] <- cop_dim_reduce(endog_cop, Z = pg_dat$Q3)
C_P_W[pg_dat$Q4==1,] <- cop_dim_reduce(endog_cop, Z = pg_dat$Q4)
# Bind to pg_dat
pg_dat <- cbind(pg_dat, C_P_W)
# Regression with copulas
ols_formula <- "YMAT ~  Price + Promo + Disp + Q2 + Q3 + Q4 + (Q1 + Q2 + Q3 + Q4):(copula1 + copula2 + copula3)"
all_endog_reg_store1 <- lm(ols_formula, data = pg_dat)
print("Price, Promotion and Display endog - store 1: Estimated coefficients")
round(coef(all_endog_reg_store1), digits = 3)

# Bootstrap results
coef_joint_test <- c("Q1:copula1","Q1:copula2","Q1:copula3","Q2:copula1","Q2:copula2","Q2:copula3","Q3:copula1","Q3:copula2","Q3:copula3","Q4:copula1","Q4:copula2","Q4:copula3")
set.seed(1)
all_endog_boot_store1 <- boot_proposed_all_endog(dat = pg_dat, B = 5000,
                                                 ols_formula = ols_formula,
                                                 coef_joint_test = coef_joint_test)
print("Price, Promotion and Display endog - store 1: Bootstrapped standarderrors")
round(all_endog_boot_store1$boot_std_error, digits = 3)
print("Price, Promotion and Display endog - store 1: P-values using bootstrapped standarderrors")
round(all_endog_boot_store1$boot_p_val_from_boot_std_error, digits = 4)
print("Price, Promotion and Display endog - store 1: Bootstrap percentile confidence intervals")
all_endog_boot_store1$boot_percentile_conf
print("Price, Promotion and Display endog - store 1: p-value from bootstrapped wald stat")
all_endog_boot_store1$wald_test_p_val
####
# Store 2
####
pg_dat <- readxl::read_excel(path = "simulated_data.xlsx", sheet = 2)
# Add quarter 1 dummy
pg_dat[,"Q1"] <- 0
pg_dat[(pg_dat$Q2 == 0 & pg_dat$Q3 == 0) & (pg_dat$Q4 == 0),"Q1"] <- 1
# List of dummy categories
categories <- list(pg_dat$Q1, pg_dat$Q2, pg_dat$Q3, pg_dat$Q4)
# Create marginal copula terms within each category
pg_dat$cop_Price <- get_copula(endog = pg_dat$Price, categories)
pg_dat$cop_Promo <- get_copula(endog = pg_dat$Promo, categories)
pg_dat$cop_Disp <- get_copula(endog = pg_dat$Disp, categories)
# Create C_P matrix
endog_cop <- list(copula1 = pg_dat$cop_Price, copula2 = pg_dat$cop_Promo, copula3 = pg_dat$cop_Disp)
C_P_W <- matrix(NA, nrow = nrow(pg_dat), ncol = length(endog_cop))
colnames(C_P_W) <- names(endog_cop)
C_P_W[pg_dat$Q1==1,] <- cop_dim_reduce(endog_cop, Z = pg_dat$Q1)
C_P_W[pg_dat$Q2==1,] <- cop_dim_reduce(endog_cop, Z = pg_dat$Q2)
C_P_W[pg_dat$Q3==1,] <- cop_dim_reduce(endog_cop, Z = pg_dat$Q3)
C_P_W[pg_dat$Q4==1,] <- cop_dim_reduce(endog_cop, Z = pg_dat$Q4)
# Bind to pg_dat
pg_dat <- cbind(pg_dat, C_P_W)
# Regression with copulas
ols_formula <- "YMAT ~  Price + Promo + Disp + Q2 + Q3 + Q4 + (Q1 + Q2 + Q3 + Q4):(copula1 + copula2 + copula3)"
all_endog_reg_store2 <- lm(ols_formula, data = pg_dat)
print("Price, Promotion and Display endog - store 2: Estimated coefficients")
round(coef(all_endog_reg_store2), digits = 3)
# Bootstrap results
coef_joint_test <- c("Q1:copula1","Q1:copula2","Q1:copula3","Q2:copula1","Q2:copula2","Q2:copula3","Q3:copula1","Q3:copula2","Q3:copula3","Q4:copula1","Q4:copula2","Q4:copula3")
set.seed(1)
all_endog_boot_store2 <- boot_proposed_all_endog(dat = pg_dat, B = 5000,
                                                 ols_formula = ols_formula,
                                                 coef_joint_test = coef_joint_test)
print("Price, Promotion and Display endog - store 2: Bootstrapped standarderrors")
round(all_endog_boot_store2$boot_std_error, digits = 3)
print("Price, Promotion and Display endog - store 2: P-values using bootstrapped standarderrors")
round(all_endog_boot_store2$boot_p_val_from_boot_std_error, digits = 4)
print("Price, Promotion and Display endog - store 2: Bootstrap percentile confidence intervals")
all_endog_boot_store2$boot_percentile_conf
print("Price, Promotion and Display endog - store 2: p-value from bootstrapped wald stat")
all_endog_boot_store2$wald_test_p_val