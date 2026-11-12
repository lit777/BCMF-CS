n = 250 # num. of obs
mcmc = 500 # num. of posterior samples
num_sim = 200 # num. of replicates


load("truth.RData")
direct_TRUE <- D_true
indirect_TRUE <- ID_true
confounder_TRUE <- 1:5

## BCMF-CS ##########################################################
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
  load(paste0("out_bcmfcs_sim5_",num,".RData"))
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
  cover_direct_bcf[,num] <- ifelse(direct_Q[1,] < direct_TRUE & direct_Q[2,] > direct_TRUE, 1, 0)
  indirect_Q <- apply(rcpp$predicted_d*rcpp$predicted_tau, 1, function(x) quantile(x, c(0.025, 0.975)))
  cover_indirect_bcf[,num] <- ifelse(indirect_Q[1,] < indirect_TRUE & indirect_Q[2,] > indirect_TRUE, 1, 0)
  
  
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


## BCMF ##########################################################
direct_BART <- matrix(nrow=n, ncol=num_sim)
indirect_BART <- matrix(nrow=n, ncol=num_sim)
cover_direct_BART <- matrix(nrow=n, ncol=num_sim)
cover_indirect_BART <-matrix(nrow=n, ncol=num_sim)
ind1_BART <- ind2_BART <- ind3_BART <- NULL
TPR <- matrix(nrow=3, ncol=num_sim)
Precision <- matrix(nrow=3, ncol=num_sim)
FPR <- matrix(nrow=3, ncol=num_sim)
F1 <- matrix(nrow=3, ncol=num_sim)

for(num in 1:num_sim){
  load(paste0("out_BCMF_sim5_",num,".RData"))
  direct_BART[,num] <- rowMeans(rcpp$predicted_zeta)
  indirect_BART[,num] <- rowMeans(rcpp$predicted_d*rcpp$predicted_tau)
  
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

  ind1_BART <- cbind(ind1_BART, colMeans(rcpp$ind))
  ind2_BART <- cbind(ind2_BART, colMeans(rcpp$ind1))
  ind3_BART <- cbind(ind3_BART, colMeans(rcpp$ind2))
  
  direct_Q <- apply(rcpp$predicted_zeta, 1, function(x) quantile(x, c(0.025, 0.975)))
  cover_direct_BART[,num] <- ifelse(direct_Q[1,] < direct_TRUE & direct_Q[2,] > direct_TRUE, 1, 0)
  indirect_Q <- apply(rcpp$predicted_d*rcpp$predicted_tau, 1, function(x) quantile(x, c(0.025, 0.975)))
  cover_indirect_BART[,num] <- ifelse(indirect_Q[1,] < indirect_TRUE & indirect_Q[2,] > indirect_TRUE, 1, 0)

}


# abs.BIAS
mean(rowMeans(apply(direct_BART, 2, function(x) abs(x-direct_TRUE))))
mean(rowMeans(apply(indirect_BART, 2, function(x) abs(x-indirect_TRUE))))

# BIAS
mean(rowMeans(apply(direct_BART, 2, function(x) (x-direct_TRUE))))
mean(rowMeans(apply(indirect_BART, 2, function(x) (x-indirect_TRUE))))

# MSE
mean(rowMeans(apply(direct_BART, 2, function(x) (x-direct_TRUE)^2)))
mean(rowMeans(apply(indirect_BART, 2, function(x) (x-indirect_TRUE)^2)))

# coverage
mean(rowMeans(cover_direct_BART))
mean(rowMeans(cover_indirect_BART))

# TPR
rowMeans(TPR)
# FPR
rowMeans(FPR)
# Precision
rowMeans(Precision)
# F1
rowMeans(F1)

## DCMF ##########################################################
direct_DART <- matrix(nrow=n, ncol=num_sim)
indirect_DART <- matrix(nrow=n, ncol=num_sim)
cover_direct_DART <- matrix(nrow=n, ncol=num_sim)
cover_indirect_DART <-matrix(nrow=n, ncol=num_sim)
ind1_DART <- ind2_DART <- ind3_DART <- NULL
TPR <- matrix(nrow=3, ncol=num_sim)
Precision <- matrix(nrow=3, ncol=num_sim)
FPR <- matrix(nrow=3, ncol=num_sim)
F1 <- matrix(nrow=3, ncol=num_sim)


for(num in 1:num_sim){
  load(paste0("out_DCMF_sim5_",num,".RData"))
  direct_DART[,num] <- rowMeans(rcpp$predicted_zeta)
  indirect_DART[,num] <- rowMeans(rcpp$predicted_d*rcpp$predicted_tau)
  
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

  
  ind1_DART <- cbind(ind1_DART, colMeans(rcpp$ind))
  ind2_DART <- cbind(ind2_DART, colMeans(rcpp$ind1))
  ind3_DART <- cbind(ind3_DART, colMeans(rcpp$ind2))
  
  direct_Q <- apply(rcpp$predicted_zeta, 1, function(x) quantile(x, c(0.025, 0.975)))
  cover_direct_DART[,num] <- ifelse(direct_Q[1,] < direct_TRUE & direct_Q[2,] > direct_TRUE, 1, 0)
  indirect_Q <- apply(rcpp$predicted_d*rcpp$predicted_tau, 1, function(x) quantile(x, c(0.025, 0.975)))
  cover_indirect_DART[,num] <- ifelse(indirect_Q[1,] < indirect_TRUE & indirect_Q[2,] > indirect_TRUE, 1, 0)
  
}


# abs.BIAS
mean(rowMeans(apply(direct_DART, 2, function(x) abs(x-direct_TRUE))))
mean(rowMeans(apply(indirect_DART, 2, function(x) abs(x-indirect_TRUE))))

# BIAS
mean(rowMeans(apply(direct_DART, 2, function(x) (x-direct_TRUE))))
mean(rowMeans(apply(indirect_DART, 2, function(x) (x-indirect_TRUE))))

# MSE
mean(rowMeans(apply(direct_DART, 2, function(x) (x-direct_TRUE)^2)))
mean(rowMeans(apply(indirect_DART, 2, function(x) (x-indirect_TRUE)^2)))

# coverage
mean(rowMeans(cover_direct_DART))
mean(rowMeans(cover_indirect_DART))

# TPR
rowMeans(TPR)
# FPR
rowMeans(FPR)
# Precision
rowMeans(Precision)
# F1
rowMeans(F1)


## LSEM ##########################################################
direct_LSEM <- matrix(nrow=n, ncol=num_sim)
indirect_LSEM <- matrix(nrow=n, ncol=num_sim)
cover_direct_LSEM <- matrix(nrow=n, ncol=num_sim)
cover_indirect_LSEM <-matrix(nrow=n, ncol=num_sim)
#nd1_LSEM <- ind2_LSEM <- ind3_LSEM <- NULL


for(num in 1:num_sim){
  load(paste0("out_LSEM_sim5_",num,".RData"))
  direct_LSEM[,num] <- out_lsem$zeta
  indirect_LSEM[,num] <- out_lsem$delta
#  ind1_DART <- cbind(ind1_DART, colMeans(rcpp$ind))
#  ind2_DART <- cbind(ind2_DART, colMeans(rcpp$ind1))
#  ind3_DART <- cbind(ind3_DART, colMeans(rcpp$ind2))
  
  direct_Q <- apply(zeta_boot_ll, 2, function(x) quantile(x, c(0.025, 0.975)))
  cover_direct_LSEM[,num] <- ifelse(direct_Q[1,] < direct_TRUE & direct_Q[2,] > direct_TRUE, 1, 0)
  indirect_Q <- apply(delta_boot_ll, 2, function(x) quantile(x, c(0.025, 0.975)))
  cover_indirect_LSEM[,num] <- ifelse(indirect_Q[1,] < indirect_TRUE & indirect_Q[2,] > indirect_TRUE, 1, 0)

}


# abs.BIAS
mean(rowMeans(apply(direct_LSEM, 2, function(x) abs(x-direct_TRUE))))
mean(rowMeans(apply(indirect_LSEM, 2, function(x) abs(x-indirect_TRUE))))

# BIAS
mean(rowMeans(apply(direct_LSEM, 2, function(x) (x-direct_TRUE))))
mean(rowMeans(apply(indirect_LSEM, 2, function(x) (x-indirect_TRUE))))

# MSE
mean(rowMeans(apply(direct_LSEM, 2, function(x) (x-direct_TRUE)^2)))
mean(rowMeans(apply(indirect_LSEM, 2, function(x) (x-indirect_TRUE)^2)))

# coverage
mean(rowMeans(cover_direct_LSEM))
mean(rowMeans(cover_indirect_LSEM))

pdf("pip_s5.pdf", height=8, width=10)
par(mfrow=c(3, 1), mai=c(0.6, 0.8, 0.3, 0.5))
plot(rowMeans(ind1_bcf), xlim=c(1, 100), ylim=c(0,1), col=c(rep("red", 5), rep("grey", 95)), ylab="PIP", xlab=expression(X[1]:X[100]))
points(rowMeans(ind1_BART), pch=2,  col=c(rep("blue", 5), rep("grey", 95)))
points(rowMeans(ind1_DART), pch=3,  col=c(rep("green", 5), rep("grey", 95)))
legend(90, 1, legend=c("BCMF-CS", "BCMF", "DCMF"), pch=c(1, 2, 3), title="Exposure Model", col=c("red", "blue", "green"))

plot(rowMeans(ind2_bcf[-1, ]), xlim=c(1, 100), ylim=c(0,1), col=c(rep("red", 5), rep("grey", 95)), ylab="PIP", xlab=expression(X[1]:X[100]))
points(rowMeans(ind2_BART[-1, ]), pch=2,  col=c(rep("blue", 5), rep("grey", 95)))
points(rowMeans(ind2_DART[-1, ]), pch=3,  col=c(rep("green", 5), rep("grey", 95)))
legend(90, 1, legend=c("BCMF-CS", "BCMF", "DCMF"), pch=c(1, 2, 3), title="Mediator Model", col=c("red", "blue", "green"))

plot(rowMeans(ind3_bcf[-(1:3), ]), xlim=c(1, 100), ylim=c(0,1), col=c(rep("red", 5), rep("grey", 95)), ylab="PIP", xlab=expression(X[1]:X[100]))
points(rowMeans(ind3_BART[-(1:3), ]), pch=2,  col=c(rep("blue", 5), rep("grey", 95)))
points(rowMeans(ind3_DART[-(1:3), ]), pch=3,  col=c(rep("green", 5), rep("grey", 95)))
legend(90, 1, legend=c("BCMF-CS", "BCMF", "DCMF"), pch=c(1, 2, 3), title="Outcome Model", col=c("red", "blue", "green"))
dev.off()
