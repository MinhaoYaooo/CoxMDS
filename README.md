# CoxMDS
`CoxMDS` is an R package to peform high dimensional mediation analysis with a survival response using Cox-proportional hazards model.  `CoxMDS` provides a finite-sample FDR control using the multiple data splitting strategy.

# Installation

You can install the development version of `CoxMDS` from Github via the `devtools` package.
```
devtools::install_github("MinhaoYaooo/CoxMDS")
```

# Detailed Tutorial

This section provides a thorough explanation of the required input data and the output objects.

## Input Data Format

The `CoxMDS` function  requires four main arguments:

* `X`: A numeric vector of length `n` (number of samples) representing the exposure variable of interest.
* `Y`: A `Surv` object for the time-to-event outcome of length `n`. This is typically created using the `survival::Surv()` function, which takes two arguments: `time` (the observed follow-up time) and `status` (the event indicator, usually 1 for event and 0 for censoring).
* `M`: A numeric matrix of size `n x p` (samples x mediators) representing the high-dimensional mediators. Each column should correspond to one mediator. It is good practice to name these columns (e.g., `M1, M2, ...`).
* `COV`: A numeric matrix of size `n x k` (samples x covariates) representing any additional covariates to be adjusted for in the model (e.g., age, sex, clinical variables). If there are no covariates, this can be omitted.

## Output Results Explanation
The CoxMDS function returns a list containing three main components:

* `$first.step`: A character vector containing the names of mediators that passed the initial screening step. This step selects mediators significantly associated with the exposure variable `X` after adjusting for covariates. These are the candidate mediators that proceed to the final selection.

* `$second.step`: A character vector containing the names of mediators that were finally selected by `CoxMDS`. These mediators have non-zero indirect effects ($\alpha \beta \neq 0$) and are considered significant discoveries with controlled FDR.

* `$estimates`: A data frame providing the estimated effect sizes for each mediator in `$second.step`. It has three columns:
  + `alpha`: The estimated coefficient for the path Exposure (`X`) -> Mediator (`M`).
  + `beta`: The estimated coefficient for the path Mediator (`M`) -> Outcome (`Y`), after adjusting for the exposure and other covariates.
  + `alpha.beta`: The product `alpha * beta`, which represents the estimated indirect effect of the exposure on the outcome through the mediator.

# Workflow Examples

The following example walks through the entire process, from data simulation to interpreting the results.

## Simulate Data

First, we load the necessary package and simulate a dataset that mimics real-world data structure.

```
# Load the package
library(CoxMDS)

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

## Run CoxMDS

Now, we apply the `CoxMDS` function to the simulated data to identify significant mediators.

```
# Perform Mediation Analysis
results <- CoxMDS(X = X, Y = Y, M = M, COV = COV, penalty = "MCP")

# View the Results
print(results)
```

## Interpret the Results

Let's break down the output list from the `results` object.


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

* Mediators listed in `$second.step` are interpreted as having statistically significant indirect effects from the exposure `X` to the survival outcome `Y` under the `CoxMDS` procedure with finite-sample FDR control.

* The `alpha.beta` value quantifies the mediator-specific indirect effect:
  + A positive `alpha.beta` suggests that mediator contributes to increased hazard (risk),
  + A negative `alpha.beta` suggests that mediator contributes to decreased hazard (protective).
