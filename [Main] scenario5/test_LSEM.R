## SIM 7

library(mediation); library(tidyverse)
rm(list=ls())

num_trial = 200
n=250; P=500; m=200; n.iter=10000

for(test_case in 1:200) {
  cat("Testing ", test_case, "of", num_trial, "at", format(Sys.time(), "%H:%M:%S"), "\n")
  # sep --------------
  
  source("source/data_simple.R")
  
  
  lsem_mediate <- function(fit_m, fit_y, formula_X,
                           data_test, mediator_name, treat_name) {
    
    X <- model.matrix(formula_X, data = data_test)
    
    gamma_m <- coef(fit_m)[str_detect(names(coef(fit_m)), treat_name)]
    gamma_y <- coef(fit_y)[str_detect(names(coef(fit_y)), treat_name)]
    xi_y <- coef(fit_y)[str_detect(names(coef(fit_y)), mediator_name)]
    
    tau_hat <- as.numeric(X %*% gamma_m)
    zeta_hat <- as.numeric(X %*% gamma_y)
    d_hat <- as.numeric(X %*% xi_y)
    delta_hat <- tau_hat * d_hat
    
    sigma_m <- sigma(fit_m)
    sigma_y <- sigma(fit_y)
    
    return(list(zeta = zeta_hat,
                delta = delta_hat,
                tau_hat = tau_hat,
                sigma_m = sigma_m,
                sigma_y = sigma_y))
  }
  
  lsem_boot <- function(fit_m, fit_y, formula_X, data_train=NULL, data_test=NULL,
                        mediator_name, treat_name, outcome_name) {
    
    formula_m <- formula(fit_m)
    formula_y <- formula(fit_y)
    
    m_fitted <- predict(fit_m)
    m_resid <- sample(resid(fit_m), replace = TRUE)
    m_boot <- m_fitted + m_resid
    data_train_boot_m <- data_train
    data_train_boot_m[[mediator_name]] <- m_boot
    fit_m_boot <- lm(formula_m, data = data_train_boot_m)
    
    y_fitted <- predict(fit_y)
    y_resid <- sample(resid(fit_y), replace = TRUE)
    y_boot <- y_fitted + y_resid
    data_train_boot_y <- data_train
    data_train_boot_y[[outcome_name]] <- y_boot
    fit_y_boot <- lm(formula_y, data = data_train_boot_y)
    
    out_boot <- lsem_mediate(fit_m_boot, fit_y_boot, formula_X,
                             data_test, mediator_name, treat_name)
    
    return(list(zeta = out_boot$zeta,
                delta = out_boot$delta))
  }
  
  train_data <- data.frame(Y_out, M_out, Y_trt, Xpred)
  
  formula_X_lsem <- as.formula(paste0("Y_trt ~ ", paste0(names(train_data)[c(4,5,6,7,8,9,10,11,12,13)], collapse = " + ")))
  formula_m_lsem <- as.formula(paste0("M_out ~ ", paste0("Y_trt*",names(train_data)[c(4,5,6,7,8,9,10,11,12,13)], collapse = " + ")))
  formula_y_lsem <- as.formula(paste0("Y_out ~ ", paste0("Y_trt*",names(train_data)[c(4,5,6,7,8,9,10,11,12,13)],"+M_out*",names(train_data)[c(4,5,6,7,8,9,10,11,12,13)], collapse = " + ")))
  fit_m_lsem <- lm(formula_m_lsem, data=train_data)
  fit_y_lsem <- lm(formula_y_lsem, data=train_data)
  
  out_lsem <- lsem_mediate(
    fit_m = fit_m_lsem, 
    fit_y = fit_y_lsem, 
    formula_X = formula_X_lsem,
    data_test = train_data,
    mediator_name = "M_out", 
    treat_name = "Y_trt"
  )
  
  zeta_hat_lsem <- out_lsem$zeta # direct
  delta_hat_lsem <- out_lsem$delta # indirect 
  sigma_m_hat_lsem <- out_lsem$sigma_m
  sigma_y_hat_lsem <- out_lsem$sigma_y
  
  out_boot_ll  <- replicate(1000, lsem_boot(fit_m_lsem, fit_y_lsem,
                                            formula_X_lsem,
                                            train_data, train_data,
                                            "M_out", "Y_trt",
                                            "Y_out"))
  zeta_boot_ll  <- matrix(unlist(out_boot_ll['zeta',]), ncol = nrow(train_data), byrow = TRUE)
  delta_boot_ll <- matrix(unlist(out_boot_ll['delta',]), ncol = nrow(train_data), byrow = TRUE)

  save(out_lsem, zeta_boot_ll,delta_boot_ll,Y_11, Y_10, Y_00, M_1, M_0, Y_trt, file=paste("out_LSEM_sim7_", test_case,".RData", sep=""))
}
