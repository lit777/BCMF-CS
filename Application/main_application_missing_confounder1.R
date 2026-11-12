#------ required libraries
library(Rcpp)
library(MCMCpack)
library(rootSolve)
library(SoftBart)
library(speff2trial)

#------ load MCMC c++ code
rm(list=ls())
sourceCpp("src/MCMC_main.cpp", rebuild = T)

#------ load the dataset
data(ACTG175)

#------ Set Treatment (TRT), Outcome (Y), Mediators (M) and Covariates (X)
Data <- ACTG175[!is.na(ACTG175$cd496), ]

Y <- Data$cd496
Trt <- Data$treat
M <- Data$cd420
X <- with(Data, cbind(age,
                      wtkg,
                      hemo,
                      homo,
                      drugs,
                      karnof,
                      oprior,
                      z30,
                      race,
                      gender,
                      str2,
                      strat,
                      symptom,
                      cd40)
)
col <- X$race+1

#------ PS estimation
P <- dim(X)[2] #<--------- Num. of Covariates
n <- dim(X)[1] #<--------- Num. of Observations

PS.fit <- glm(Trt~X, family=binomial()) # to prevent perfect separation
PS <- predict(PS.fit, type="response")

#------ (M1,M0) estimation
M.fit <- lm(M~Trt+X)
Mest1 <- predict(M.fit, newdata = data.frame(Trt=rep(1,n)))
Mest0 <- predict(M.fit, newdata = data.frame(Trt=rep(0,n)))

X.ps <- quantile_normalize_bart(cbind(PS, X))
X.Mps <- quantile_normalize_bart(cbind(Mest1, Mest0,PS, X))
X.M <- quantile_normalize_bart(cbind(Mest1, Mest0, X))

x.index <- sample(c(0, 100, 200, 300, 400, 500, 600, 700, 800), n, replace=T)
X.M_mult <- quantile_normalize_bart(rbind(cbind(Mest1, Mest0, X[,1:14], x.index)))
X_mult <- quantile_normalize_bart(rbind(cbind(PS, X[,1:14], x.index)))

X <- quantile_normalize_bart(X)



#------ MCMC settings
n.iter=30000;

nu <- 3    # default setting (nu, q) = (3, 0.90) from Chipman et al. 2010
m <- 100                  # Num. of trees
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


#------ Main MCMC running
rcpp = MCMC(X, X.ps,  X.M, X.Mps, X_mult,X.M_mult, Trt, M, Y,
            p.grow, p.prune, p.change, m, m, 20, m, 20, 20, nu,
            lambda, lambda_m, lambda_y, 0.1,0.1,0.1,0.1, alpha, alpha_modifier, beta, n.iter)


#------ Plots
par(mfrow=c(3,1))
plot(rowMeans(rcpp$predicted_Y11-rcpp$predicted_Y00), col=col, pch=col, main="Total Effects", ylab=expression(Y["11"]-Y["00"]), xlab="Patients ID")
plot(rowMeans(rcpp$predicted_Y10-rcpp$predicted_Y00), col=col, pch=col, main="Direct Effects", ylab=expression(Y["10"]-Y["00"]), xlab="Patients ID")
plot(rowMeans(rcpp$predicted_Y11-rcpp$predicted_Y10), col=col, pch=col, main="Indirect Effects", ylab=expression(Y["11"]-Y["10"]), xlab="Patients ID")
dev.off()

var_names <- c("age", "wtkg",
               "hemo", "homo",
               "drugs", "karnof",
               "oprior", "z30",
               "race",
               "gender", "str2",
               "strat", "symptom", "cd40")
par(mar = c(7, 4, 4, 2) + 0.1)

plot(colMeans(rcpp$ind),
     ylim = c(0, 1),
     xaxt = "n",        
     yaxt = "n",         
     xlab = "",
     ylab = "PIP",
     pch  = 1)

points(colMeans(rcpp$ind1)[-1],         pch = 2)
points(colMeans(rcpp$ind2)[-c(1,2,3)],  pch = 3)

axis(2, at = c(0, 0.5, 1), las = 1)
abline(h = 0.5, lty = 2, col = "gray60")  

axis(1, at = seq_along(var_names), labels = FALSE)


text(x      = seq_along(var_names),
     y      = par("usr")[3] - 0.03,  
     labels = var_names,
     srt    = 45,       
     adj    = c(1, 1),  
     xpd    = TRUE,      
     cex    = 0.8)

legend("topright",
       legend = c("Exposure", "Mediator", "Outcome"),
       pch    = c(1, 2, 3),
       bty    = "n")
