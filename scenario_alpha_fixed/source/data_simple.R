
cov <- list()
for(i in 1:2){
  cov[[i]] <- runif(n,0,1)
}
cov[[3]] <- rnorm(n, 0, 1)
cov[[4]] <- rnorm(n, 0, 1)
cov[[5]] <- rnorm(n, 0, 1) 
for(i in 6:P){
  cov[[i]] <- rnorm(n,0,1)
}

Xpred <- do.call(cbind, cov)
Xcut <- lapply(1:dim(Xpred)[2], function(t) sort(unique(Xpred[,t]))) # e.g. unique values of predictors


MU_func <- -0.3+Xpred[,1]-1.5*Xpred[,2]+0.8*abs(Xpred[,3]-1)+0.5*Xpred[,4]*Xpred[,5]

Y_trt <- rbinom(n, 1, exp(MU_func)/(1+exp(MU_func))) # exposure model

M_err <- rnorm(n,0,0.5)
M_1 <- 1+0.9*1+1*(Xpred[,1])-0.5*Xpred[,2]+0.5*abs(Xpred[,3]-1) - 0.5*Xpred[,4]^2+sin(Xpred[,5])+M_err
M_0 <- 1+0.9*0+1*(Xpred[,1])-0.5*Xpred[,2]+0.5*abs(Xpred[,3]-1) - 0.5*Xpred[,4]^2+sin(Xpred[,5])+M_err
M_out <- ifelse(Y_trt==1, M_1, M_0)



Y_err <- rnorm(n, 0, 0.3)
Y_11 <- 1*(Xpred[,1])+1.5*Xpred[,2]+0.5*Xpred[,3]*Xpred[,4] - 0.4*Xpred[,5]^2-1.5*(M_1)+1*1+Y_err
Y_00 <- 1*(Xpred[,1])+1.5*Xpred[,2]+0.5*Xpred[,3]*Xpred[,4] - 0.4*Xpred[,5]^2-1.5*(M_0)+1*0+Y_err
Y_10 <- 1*(Xpred[,1])+1.5*Xpred[,2]+0.5*Xpred[,3]*Xpred[,4] - 0.4*Xpred[,5]^2-1.5*(M_0)+1*1+Y_err

Y_out <- ifelse(Y_trt==1, Y_11, Y_00)


### Initial Setup (priors, initial values and hyper-parameters)
p.grow <- 0.28            # Prob. of GROW
p.prune <- 0.28           # Prob. of PRUNE
p.change <- 0.44          # Prob. of CHANGE
# m <- 100                  # Num. of Trees: default setting 100
# sigma2_m <- var(M_out)       # Initial value of SD^2
# sigma2_y <- var(Y_out) 
sigma2_m <- 1
sigma2_y <- 1 


nu <- 3                   # default setting (nu, q) = (3, 0.90) from Chipman et al. 2010
f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2_y)
lambda_y <- rootSolve::uniroot.all(f, c(0.1^5,10))
f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2_m)
lambda_m <- rootSolve::uniroot.all(f, c(0.1^5,10))

sigma2 <- 1
f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2)
lambda <- rootSolve::uniroot.all(f, c(0.1^5,10))

alpha <- 0.95             # alpha (1+depth)^{-beta} where depth=0,1,2,...
alpha_modifier <- 0.5             # alpha (1+depth)^{-beta} where depth=0,1,2,...
beta <- 2                 # default setting (alpha, beta) = (0.95, 2)


