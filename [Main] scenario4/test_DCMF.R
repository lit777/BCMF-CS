## SIM 6

library(Rcpp); library(MCMCpack); library(rootSolve); library(SoftBart)
rm(list=ls())
sourceCpp("src/MCMC_main_DCMF.cpp", rebuild = T)

num_trial = 200
n=250; P=100; m=200; n.iter=10000



for(test_case in 1:200) {
  cat("Testing ", test_case, "of", num_trial, "at", format(Sys.time(), "%H:%M:%S"), "\n")
  
  source("source/data_simple.R")
  
  
    PS.fit <- glm(Y_trt~Xpred[,1:10], family=binomial()) # to prevent perfect separation
    PS <- predict(PS.fit, type="response")
    M.fit <- lm(M_out~Y_trt+Xpred)
    Mest1 <- predict(M.fit, newdata = data.frame(Y_trt=rep(1,n)))
    Mest0 <- predict(M.fit, newdata = data.frame(Y_trt=rep(0,n)))
    
    Xpred.ps <- cbind(PS, quantile_normalize_bart(Xpred))
    Xpred.Mps <- cbind(quantile_normalize_bart(cbind(Mest1, Mest0)),PS, quantile_normalize_bart(Xpred))
    Xpred.M <- quantile_normalize_bart(cbind(Mest1, Mest0, Xpred))
    Xpred <- quantile_normalize_bart(Xpred)
    
    Xpred.M_mult <- Xpred.M
    Xpred_mult <- Xpred
   
  rcpp = MCMC(Xpred, Xpred.ps,  Xpred.M, Xpred.Mps, Xpred_mult, Xpred.M_mult, Y_trt, M_out, Y_out, p.grow, p.prune, p.change, m, m, 20, m, 20, 20, nu, lambda, lambda_m, lambda_y, 0.5,0.5,0.5,0.5,0.5,0.5, alpha, alpha_modifier, beta, n.iter)

  save(rcpp, Y_11, Y_10, Y_00, M_1, M_0, Y_trt, file=paste("out_DCMF_sim6_", test_case,".RData", sep=""))
}
