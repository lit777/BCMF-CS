## SIM 2

library(Rcpp); library(MCMCpack); library(rootSolve); library(SoftBart)
rm(list=ls())
sourceCpp("src/MCMC_main_test_binary_fixed.cpp", rebuild = T)

num_trial = 200
n=250; P=100; m=200; n.iter=10000



for(test_case in 1:100) {
  cat("Testing ", test_case, "of", num_trial, "at", format(Sys.time(), "%H:%M:%S"), "\n")
  # sep --------------
  
  source("source/data_simple.R")
  
 
  
#  f <- function(scale) qcauchy(0.75, 0, scale) - 2*sd(M_out-mean(M_out))   # first
  f <- function(scale) qcauchy(0.75, 0, scale) - 2*1   # first
  sigma_mu_m_mu_sigma <- uniroot.all(f, c(0.1^5, 100))
  
#  f <- function(sd) qnorm(0.75, 0, sd) - sd(M_out-mean(M_out))             # second
  f <- function(sd) qnorm(0.75, 0, sd) - 1            # second
  sigma_mu_m_tau_sigma <- uniroot.all(f, c(0.1^5, 100))
  
#  f <- function(scale) qcauchy(0.75, 0, scale) - 2*sd(Y_out-mean(Y_out))   # first
  f <- function(scale) qcauchy(0.75, 0, scale) - 2*1   # first
  sigma_mu_y_mu_sigma <- uniroot.all(f, c(0.1^5, 100))
  
#  f <- function(sd) qnorm(0.75, 0, sd) - sd(Y_out-mean(Y_out))             # second
  f <- function(sd) qnorm(0.75, 0, sd) - 1             # second
  sigma_mu_y_tau_sigma <- uniroot.all(f, c(0.1^5, 100))
  
    PS.fit <- glm(Y_trt~Xpred[,1:10], family=binomial()) # to prevent perfect separation
    PS <- predict(PS.fit, type="response")
    M.fit <- lm(M_out~Y_trt+Xpred)
    Mest1 <- predict(M.fit, newdata = data.frame(Y_trt=rep(1,n)))
    Mest0 <- predict(M.fit, newdata = data.frame(Y_trt=rep(0,n)))
    
    Xpred.ps <- cbind(PS, quantile_normalize_bart(Xpred))
    Xpred.Mps <- cbind(quantile_normalize_bart(cbind(Mest1, Mest0)),PS, quantile_normalize_bart(Xpred))
    Xpred.M <- quantile_normalize_bart(cbind(Mest1, Mest0, Xpred))
    Xpred <- quantile_normalize_bart(Xpred)
    
    Xpred_mult <- Xpred.M
    # 
  rcpp = MCMC(Xpred, Xpred.ps,  Xpred.M, Xpred.Mps, Xpred_mult, Y_trt, M_out, Y_out, p.grow, p.prune, p.change, m, m, 20, m, 20, 20, nu, lambda, lambda_m, lambda_y, 0.5,0.5,0.5,0.5, alpha, alpha_modifier, beta, n.iter)

  save(rcpp, Y_11, Y_10, Y_00, M_1, M_0, Y_trt, file=paste("out_bcf_alpha_fixed_", test_case,".RData", sep=""))
}
