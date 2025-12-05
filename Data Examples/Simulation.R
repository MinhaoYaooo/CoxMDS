library(CoxMDS)
library(CoxMKF)
library(MASS)


N <- 500
n <- 500                                   #number of replications
p <- 2000                                  #number of mediators
theta1 <- 0.3                              #coefficients(Z1-->M)
theta2 <- 0.2                              #coefficients(Z2-->M)
phi <- c(0.3, -0.2)                        #coefficients(covariates-->outcome)
rho <- 0.8 # correlation
c0 <- 10 # control for censoring rate
alpha=rep(0,p)                             #coefficients (mediator~exposure)
beta=rep(0,p)                              #coefficients (outcome~mediators)
alpha[1:12] <- 0.50*c(1,1,1,1,-1,-1,-1,-1,1,-1,0,0)
beta[1:12] <- 0.35*c(1,1,1,1,1,1,1,1,0,0,1,-1)


sim_data <- function(n, p, c0){
  X <- t(t(rbinom(n, 1, 0.6)))               #exposure
  Z1 <- t(t(rbinom(n, 1, 0.3)))              #covariates Z1
  Z2 <- t(t(runif(n, 0, 1)))                 #covariates Z2
  Z <- cbind(Z1, Z2)
  ck <- t(runif(p, 0, 1))
  
  Sigma <- matrix(nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i,j] <- rho^abs(i-j)
    }
  }
  
  e <- mvrnorm(n=n, mu=rep(0,p), Sigma = Sigma)
  M <- matrix(rep(ck, n),ncol = p, byrow = TRUE) + X%*%alpha + Z[,1]*theta1+Z[,2]*theta2+ e
  colnames(M) <- paste0("M", 1:ncol(M))
  haz <- 0.5*exp(0.5*X+0.3*Z[,1]-0.2*Z[,2]+M%*%beta)   #baseline hazard function lambda0 <- 0.5
  ft <- rexp(n, haz)
  ct <- rexp(n, c0)                         #censoring time
  time <- pmin(ft, ct)                       #observed time
  status <- as.numeric(ft <= ct)             #censoring indicator
  Y <- survival::Surv(time, status)
  COV <- Z
  return(list(Y=Y, X=X, M=M, COV=COV, status=status))
}

true_set <- colnames(sim.dat$M)[which(alpha!=0 & beta!=0)]
false_set <- colnames(sim.dat$M)[-which(alpha!=0 & beta!=0)]

FDR.MDS <- rep(NA, N); TPP.MDS <- rep(NA,N)
FDR.MKF <- rep(NA, N); TPP.MKF <- rep(NA,N)

for (i in 1:N) {
  sim.dat <- sim_data(n=n, p=p, c0=c0)
  
  results.CoxMDS <- CoxMDS(sim.dat$X, sim.dat$Y, sim.dat$M, sim.dat$COV, penalty = 'MCP', q2=0.1)
  results.CoxMKF <- CoxMKF(sim.dat$X, sim.dat$Y, sim.dat$M, sim.dat$COV, penalty = 'MCP', q2=0.1)
  
  ID.MDS <- results.CoxMDS$second.step
  FDR.MDS[i] <- length(intersect(false_set, ID.MDS)) / length(ID.MDS)
  TPP.MDS[i] <- length(intersect(true_set, ID.MDS)) / length(true_set)
  
  ID.MKF <- results.CoxMKF$second.step
  FDR.MKF[i] <- length(intersect(false_set, ID.MKF)) / length(ID.MKF)
  TPP.MKF[i] <- length(intersect(true_set, ID.MKF)) / length(true_set)
}
  


