n = 250 # num. of obs
mcmc = 500 # num. of posterior samples
num_sim = 200 # num. of replicates


load("truth.RData")
direct_TRUE <- D_true
indirect_TRUE <- ID_true

## BCMF-CS (50) ##########################################################
direct_bcf <- matrix(nrow=n, ncol=num_sim)
indirect_bcf <-matrix(nrow=n, ncol=num_sim)
cover_direct_bcf <- matrix(nrow=n, ncol=num_sim)
cover_indirect_bcf <-matrix(nrow=n, ncol=num_sim)
ind1_bcf <- ind2_bcf <- ind3_bcf <- NULL

for(num in 1:num_sim){
  load(paste0("out_bcmfcs_sim_",num,".RData"))
  direct_bcf[,num] <- rowMeans(rcpp$predicted_zeta)
  indirect_bcf[,num] <- rowMeans(rcpp$predicted_d*rcpp$predicted_tau)
  ind1_bcf <- cbind(ind1_bcf, colMeans(rcpp$ind))
  ind2_bcf <- cbind(ind2_bcf, colMeans(rcpp$ind1))
  ind3_bcf <- cbind(ind3_bcf, colMeans(rcpp$ind2))
  
  direct_Q <- apply(rcpp$predicted_zeta, 1, function(x) quantile(x, c(0.025, 0.975)))
  cover_direct_bcf[,num] <- ifelse(direct_Q[1,] < -1 & direct_Q[2,] > -1, 1, 0) 
  indirect_Q <- apply(rcpp$predicted_d*rcpp$predicted_tau, 1, function(x) quantile(x, c(0.025, 0.975)))
  cover_indirect_bcf[,num] <- ifelse(indirect_Q[1,] < -1.95 & indirect_Q[2,] > -1.95, 1, 0) 
  
}

mean(abs(direct_bcf-direct_TRUE)) # BIAS
mean((direct_bcf-direct_TRUE)^2) # MSE
mean(abs(indirect_bcf-indirect_TRUE)) # BIAS
mean((indirect_bcf-indirect_TRUE)^2) # MSE
(rowMeans(ind1_bcf))
(rowMeans(ind2_bcf))
(rowMeans(ind3_bcf))
rowMeans(cover_direct_bcf)
rowMeans(cover_indirect_bcf)

## BCMF-CS (100) ##########################################################
direct_bcf <- matrix(nrow=n, ncol=num_sim)
indirect_bcf <-matrix(nrow=n, ncol=num_sim)
cover_direct_bcf <- matrix(nrow=n, ncol=num_sim)
cover_indirect_bcf <-matrix(nrow=n, ncol=num_sim)
ind1_bcf <- ind2_bcf <- ind3_bcf <- NULL

for(num in 1:num_sim){
  load(paste0("out_bcmfcs_sim_100_",num,".RData"))
  direct_bcf[,num] <- rowMeans(rcpp$predicted_zeta)
  indirect_bcf[,num] <- rowMeans(rcpp$predicted_d*rcpp$predicted_tau)
  ind1_bcf <- cbind(ind1_bcf, colMeans(rcpp$ind))
  ind2_bcf <- cbind(ind2_bcf, colMeans(rcpp$ind1))
  ind3_bcf <- cbind(ind3_bcf, colMeans(rcpp$ind2))
  
  direct_Q <- apply(rcpp$predicted_zeta, 1, function(x) quantile(x, c(0.025, 0.975)))
  cover_direct_bcf[,num] <- ifelse(direct_Q[1,] < -1 & direct_Q[2,] > -1, 1, 0) 
  indirect_Q <- apply(rcpp$predicted_d*rcpp$predicted_tau, 1, function(x) quantile(x, c(0.025, 0.975)))
  cover_indirect_bcf[,num] <- ifelse(indirect_Q[1,] < -1.95 & indirect_Q[2,] > -1.95, 1, 0) 
  
}

mean(abs(direct_bcf-direct_TRUE)) # BIAS
mean((direct_bcf-direct_TRUE)^2) # MSE
mean(abs(indirect_bcf-indirect_TRUE)) # BIAS
mean((indirect_bcf-indirect_TRUE)^2) # MSE
(rowMeans(ind1_bcf))
(rowMeans(ind2_bcf))
(rowMeans(ind3_bcf))
rowMeans(cover_direct_bcf)
rowMeans(cover_indirect_bcf)

