n = 250 # num. of obs
mcmc = 500 # num. of posterior samples
num_sim = 200 # num. of replicates

# c(1,2, 3, 4, 5)
confounder_TRUE <- c(1, 2, 3, 4, 5)


direct_TRUE <- rep(1, n)
indirect_TRUE <- rep(-1.35, n)

## BCF ##########################################################
direct_bcf <- matrix(nrow=n, ncol=num_sim)
indirect_bcf <-matrix(nrow=n, ncol=num_sim)
cover_direct_bcf <- matrix(nrow=n, ncol=num_sim)
cover_indirect_bcf <-matrix(nrow=n, ncol=num_sim)
ind1_bcf <- ind2_bcf <- ind3_bcf <- NULL
TPR <- matrix(nrow=3, ncol=num_sim)
Precision <- matrix(nrow=3, ncol=num_sim)
FPR <- matrix(nrow=3, ncol=num_sim)
F1 <- matrix(nrow=3, ncol=num_sim)

for(num in 1:num_sim){
  load(paste0("out_bcf_xi_xj_",num,".RData"))
  direct_bcf[,num] <- rowMeans(rcpp$predicted_zeta)
  indirect_bcf[,num] <- rowMeans(rcpp$predicted_d*rcpp$predicted_tau)
  
  TP1 <- length(intersect(which(colMeans(rcpp$ind) > 0.5), confounder_TRUE))
  TP2 <- length(intersect(which(colMeans(rcpp$ind1[,-1]) > 0.5), confounder_TRUE))
  TP3 <- length(intersect(which(colMeans(rcpp$ind2[,-(1:3)]) > 0.5), confounder_TRUE))
  FP1 <- length(which(colMeans(rcpp$ind) > 0.5)) - TP1
  FP2 <- length(which(colMeans(rcpp$ind1[,-1]) > 0.5)) - TP2
  FP3 <- length(which(colMeans(rcpp$ind2[,-(1:3)]) > 0.5)) - TP3
  
  TPR[1, num] <- TP1 / (length(confounder_TRUE))
  TPR[2, num] <- TP2 / (length(confounder_TRUE))
  TPR[3, num] <- TP3 / (length(confounder_TRUE))
  Precision[1, num] <- TP1 / (TP1 + FP1)
  Precision[2, num] <- TP2 / (TP2 + FP2)
  Precision[3, num] <- TP3 / (TP3 + FP3)
  FPR[1, num] <- FP1 / (100-length(confounder_TRUE))
  FPR[2, num] <- FP2 / (100-length(confounder_TRUE))
  FPR[3, num] <- FP3 / (100-length(confounder_TRUE))
  F1[1, num] <- (2*Precision[1, num]*TPR[1, num]) / (Precision[1, num] + TPR[1, num])
  F1[2, num] <- (2*Precision[2, num]*TPR[2, num]) / (Precision[2, num] + TPR[2, num])
  F1[3, num] <- (2*Precision[3, num]*TPR[3, num]) / (Precision[3, num] + TPR[3, num])

  
  ind1_bcf <- cbind(ind1_bcf, colMeans(rcpp$ind))
  ind2_bcf <- cbind(ind2_bcf, colMeans(rcpp$ind1))
  ind3_bcf <- cbind(ind3_bcf, colMeans(rcpp$ind2))
  
  
  
  direct_Q <- apply(rcpp$predicted_zeta, 1, function(x) quantile(x, c(0.025, 0.975)))
  cover_direct_bcf[,num] <- ifelse(direct_Q[1,] < 1 & direct_Q[2,] > 1, 1, 0) 
  indirect_Q <- apply(rcpp$predicted_d*rcpp$predicted_tau, 1, function(x) quantile(x, c(0.025, 0.975)))
  cover_indirect_bcf[,num] <- ifelse(indirect_Q[1,] < -1.35 & indirect_Q[2,] > -1.35, 1, 0) 
  
}

# abs.BIAS
mean(rowMeans(apply(direct_bcf, 2, function(x) abs(x-direct_TRUE))))
mean(rowMeans(apply(indirect_bcf, 2, function(x) abs(x-indirect_TRUE))))


# BIAS
mean(rowMeans(apply(direct_bcf, 2, function(x) (x-direct_TRUE))))
mean(rowMeans(apply(indirect_bcf, 2, function(x) (x-indirect_TRUE))))

# MSE
mean(rowMeans(apply(direct_bcf, 2, function(x) (x-direct_TRUE)^2)))
mean(rowMeans(apply(indirect_bcf, 2, function(x) (x-indirect_TRUE)^2)))

# coverage
mean(rowMeans(cover_direct_bcf))
mean(rowMeans(cover_indirect_bcf))

# TPR
rowMeans(TPR)
# FPR
rowMeans(FPR)
# Precision
rowMeans(Precision)
# F1
rowMeans(F1)
