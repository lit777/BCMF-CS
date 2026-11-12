#------ required libraries
library(Rcpp)
library(MCMCpack)
library(rootSolve)
library(SoftBart)

#------ load MCMC c++ code
rm(list=ls())
sourceCpp("src/MCMC_main_test_binary.cpp", rebuild = T)

#------ load the dataset (for a 50km radius analysis)
load("source/Master2014_PM25_150km_annual.RData")

#------ Set Treatment (TRT), Outcome (Y), Mediators (M) and Covariates (X)
Data <- Master
Data <- subset(Data, !is.na(PM25))
Data <- subset(Data, !is.na(RHUM))
Data <- subset(Data, !is.na(Temperature))
Data <- subset(Data, !is.na(APCP))
Data <- subset(Data, totHeatInput!=0)
Data <- subset(Data, totLoad!=0) # this process automatically excludes power plants with totHeatInput==0
Data <- subset(Data, !is.na(Sulfur_Content))
Data <- subset(Data, !is.na(pctCapacity))

Y <- Data$PM25
Trt <- ifelse(Data$SO2.sc < 0.5, 0, 1)
M <- log(Data$totSO2emissions)
X <- as.matrix(cbind(with(Data, cbind(pctS_n_CR,
                      NumNOxControls,
                      Sulfur_Content,
                      totNOxemissions_pre,
                      totCO2emissions_pre,
                      totLoad,
                      totHeatInput,
                      pctCapacity,
                      Phase2,
#                      RHUM,
#                      APCP,
                      Population,
                      Income,
#                      Temperature,
                      totOpTime)), Data[,27:242]
))


#------ PS estimation

P <- dim(X)[2] #<--------- Num. of Covariates
n <- dim(X)[1] #<--------- Num. of Observations

PS.fit <- glm(Trt~X[,1:12], family=binomial()) # to prevent perfect separation
PS <- predict(PS.fit, type="response")
M.fit <- lm(M~Trt+X)
Mest1 <- predict(M.fit, newdata = data.frame(Trt=rep(1,n)))
Mest0 <- predict(M.fit, newdata = data.frame(Trt=rep(0,n)))

X.ps <- cbind(PS, quantile_normalize_bart(X))
X.Mps <- cbind(quantile_normalize_bart(cbind(Mest1, Mest0)),PS, quantile_normalize_bart(X))
X.M <- quantile_normalize_bart(cbind(Mest1, Mest0, X))
X <- quantile_normalize_bart(X)

X.M_mult <- X.M
X_mult <- X


#------ MCMC settings
n.iter=50000; 

nu <- 3    # default setting (nu, q) = (3, 0.90) from Chipman et al. 2010
m <- 200                  # Num. of trees
p.grow <- 0.28            # Prob. of GROW
p.prune <- 0.28           # Prob. of PRUNE
p.change <- 0.44          # Prob. of CHANGE

sigma2_m <- 1        # Initial value of SD^2
sigma2_y <- 1

f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2_y)
lambda_y <- rootSolve::uniroot.all(f, c(0.1^5,10))

f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2_m)
lambda_m <- rootSolve::uniroot.all(f, c(0.1^5,10))

sigma2 <- 1
f <- function(lambda) invgamma::qinvgamma(0.90, nu/2, rate = lambda*nu/2, lower.tail = TRUE, log.p = FALSE) - sqrt(sigma2)
lambda <- rootSolve::uniroot.all(f, c(0.1^5,10))


alpha <- 0.95             # alpha (1+depth)^{-beta} where depth=0,1,2,...
beta <- 2                 # default setting (alpha, beta) = (0.95, 2)
alpha_modifier <- 0.5             # alpha (1+depth)^{-beta} where depth=0,1,2,...

f <- function(scale) qcauchy(0.75, 0, scale) - 2*1
sigma_mu_m_mu_sigma <- uniroot.all(f, c(0.1^5, 100))
f <- function(sd) qnorm(0.75, 0, sd) - 1
sigma_mu_m_tau_sigma <- uniroot.all(f, c(0.1^5, 100))
f <- function(scale) qcauchy(0.75, 0, scale) - 2*1
sigma_mu_y_mu_sigma <- uniroot.all(f, c(0.1^5, 100))
f <- function(sd) qnorm(0.75, 0, sd) - 1
sigma_mu_y_tau_sigma <- uniroot.all(f, c(0.1^5, 100))


#------ Main MCMC running
rcpp = MCMC(X, X.ps,  X.M, X.Mps, X_mult,X.M_mult, Trt, M, Y, 
            p.grow, p.prune, p.change, m, m, 20, m, 20, 20, nu, 
            lambda, lambda_m, lambda_y, 0.1,0.1,0.1,0.1, alpha, alpha_modifier, beta, n.iter)




