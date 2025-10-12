# CoxMDS
`CoxMDS` is an R package to peform high dimensional mediation analysis with a survival response using Cox-proportional hazards model.  `CoxMDS` provides a finite-sample FDR control using the multiple data splitting strategy.

# Installation

You can install the development version of `CoxMDS` from Github via the `devtools` package.
```
devtools::install_github("MinhaoYaooo/CoxMDS")
```

# Examples
First, we generate the data:
```
n <- 500                                   #number of samples
p <- 1000                                  #number of mediators
alpha=rep(0,p)                             #coefficients (mediator~exposure)
beta=rep(0,p)                              #coefficients (outcome~mediators)
alpha[1:12] <- c(0.55,0.45,-0.4,-0.45,0.5,0.6,-0.4,-0.46,-0.4,0.5,0,0)
beta[1:12] <- c(0.52,0.45,0.4,0.4,-0.54,-0.6,-0.4,-0.5,0,0,0.4,-0.8)
X <- t(t(rbinom(n, 1, 0.6)))               #exposure
Z1 <- t(t(rbinom(n, 1, 0.3)))              #covariates Z1
theta1 <- 0.3                              #coefficients(Z1-->M)
Z2 <- t(t(runif(n, 0, 1)))                 #covariates Z2
theta2 <- 0.2                              #coefficients(Z2-->M)
Z <- cbind(Z1, Z2)
phi <- c(0.3, -0.2)                        #coefficients(covariates-->outcome)
ck <- t(runif(p, 0, 1))
M <- matrix(0, n, p)                       #mediators
for(i in 1:n){
  e <- rnorm(p, sd = 1)
  M[i,] <- ck+X[i]*alpha+Z[i,1]*theta1+Z[i,2]*theta2+e
}
colnames(M) <- paste0("M", 1:ncol(M))
haz <- 0.5*exp(0.5*X+0.3*Z[,1]-0.2*Z[,2]+M%*%beta)   #baseline hazard function lambda0 <- 0.5
ft <- rexp(n, haz)
ct <- rexp(n, 0.7)                         #censoring time
time <- pmin(ft, ct)                       #observed time
status <- as.numeric(ft <= ct)             #censoring indicator
Y <- survival::Surv(time, status)
COV <- Z
```
Then we apply `CoxMDS` to select the mediators:
```
results <- CoxMDS(X, Y, M, COV, penalty="MCP")
```

# Results

Here is an example of the results:

```
$first.step
 [1] "M1"   "M2"   "M3"   "M4"   "M5"   "M6"   "M7"   "M8"   "M9"  
[10] "M10"  "M17"  "M417" "M496" "M919" "M978"

$second.step
[1] "M1" "M2" "M3" "M4" "M5" "M6" "M7" "M8"

$estimates
        alpha       beta alpha.beta
M1  0.6596348  0.4083061  0.2693329
M2  0.4783746  0.4450930  0.2129212
M3 -0.3939130  0.3636737 -0.1432558
M4 -0.5280515  0.3384643 -0.1787266
M5  0.3670641 -0.4475635 -0.1642845
M6  0.5909356 -0.4930380 -0.2913537
M7 -0.5962105 -0.3871237  0.2308072
M8 -0.4476374 -0.4785410  0.2142129
```

The `$first.step` returns candidate mediators after X->M screening. The `$second.step` returns mediators selected by CoxMDS with non-zero indirect effects. The `$estimates` returns estimates of the effects.


