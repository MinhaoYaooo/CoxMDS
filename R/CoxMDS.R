#' @title Screening Based on X->M with FDR control
#' @description The first step of CoxMDS, which conducts preliminary screening with the p-values from X to M with FDR control.
#' @param X a vector of the exposure.
#' @param M a matrix of the mediators.
#' @param COV optional, a matrix of the potential covariates.
#' @param intercept binary, whether an intercept should be included in the regression M~X, the default value is TRUE.
#' @param fdr the pre-defined FDR level in the screening step, the default value is 0.2.
#'
#' @return A list of the screening results including the IDs, pvalues and coefficients.
#' @export

screen_FDR <- function(X, M, COV, intercept=TRUE, fdr=0.2){
  p <- ncol(M)
  p_alpha <- rep(NA, p)
  coef_alpha <- rep(NA, p)
  for(j in 1:p){
    if (is.null(COV)) {
      MX <- data.frame(M = M[, j], X = X)
    } else {
      MX <- data.frame(M = M[, j], X = X, COV = COV)
    }
    if(intercept){
      fit <- stats::lm(M ~1+., data = MX)
      p_alpha[j] <- summary(fit)$coef[2,4]   #p-value for alpha
      coef_alpha[j] <- summary(fit)$coef[2,1]
    }else{
      fit <- stats::lm(M ~0+., data = MX)
      p_alpha[j] <- summary(fit)$coef[1,4]   #p-value for alpha
      coef_alpha[j] <- summary(fit)$coef[1,1]
    }
  }
  names(p_alpha) <- colnames(M); names(coef_alpha) <- colnames(M)
  p_alpha.adjust <- p.adjust(p_alpha, method = 'fdr')
  ID_screen <- which(p_alpha.adjust<=fdr)
  return(list(ID_screen=ID_screen, pval=p_alpha, coef=coef_alpha))
}


#' @title Data Split for Survival Response
#' @describeIn Split the data for survival outcome and compute the corresponding statistics.
#' @param X a vector of the exposure.
#' @param Y a vector of the outcome.
#' @param MS a matrix of the mediators after the pre-liminary screening.
#' @param COV optional, a matrix of the potential covariates.
#' @param penalty the penalty used in the cox regression.
#' @param q the pre-defined FDR level in the data splitting procedure.
#' @param n_split the pre-defined number of multiple data splittings.
#'
#' @return A vector of mediator IDs selected by AKO.
#' @export

datasplit_surv <- function(X, Y, MS, COV, penalty=c('MCP','lasso','SCAD'),  q=0.1, n_split = 25){

  pS <- ncol(MS); n <- nrow(MS)
  inclusion_rate = matrix(0, n_split, pS)
  n_select <- rep(0, n_split)

  for (i in 1:n_split) {

    ### randomly split the data
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)

    ### fit penalized cox regression on the first half of the data
    XM <- cbind(MS, X)
    if (is.null(COV)) {
      fit <- ncvreg::ncvsurv(XM[sample_index1,], Y[sample_index1],
                             penalty = penalty, nlambda = 5000,
                             penalty.factor = c(rep(1, pS), 0))
    }else {
      XM_COV <- cbind(XM, COV)
      fit <-  ncvreg::ncvsurv(XM_COV[sample_index1,], Y[sample_index1],
                              penalty = penalty, nlambda = 1000,
                              penalty.factor = c(rep(1, pS), rep(0, 1 + ncol(COV))))
    }

    lam <- fit$lambda[which.min(BIC(fit))]
    Coefficients <- coef(fit, lambda = lam)
    beta1 <- Coefficients[1:pS]
    nonzero_index <- which(beta1 != 0)
    zero_index <- which(beta1 == 0)

    ### fit cox regression on the second half of the data
    if(length(nonzero_index)!=0){
      if(is.null(COV)){
        DATA = data.frame(Y=Y, XM[,-zero_index])
      }else{
        DATA = data.frame(Y=Y, XM_COV[,-zero_index])
      }
      fit2 <- survival::coxph(Y ~ ., data = DATA[sample_index2,])
      beta2 <- rep(0,pS)
      beta2[nonzero_index] <- summary(fit2)$coefficients[1:length(nonzero_index)]

      ### calculate the mirror statistics
      mirror_stat <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))

      selected_index = analys(mirror_stat, q)
      DS_selected_index = selected_index

      if(length(selected_index)!=0){
        n_select[i] <- length(selected_index)
        inclusion_rate[i, selected_index] <- 1/n_select[i]
      }

    }else{
      DS_selected_index = NULL
    }
  }


  ### aggregate the multiple data splitting results
  inclusion_rate <- apply(inclusion_rate, 2, mean)

  feature_rank <- order(inclusion_rate)
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  if(length(feature_rank)!=0){
    null_feature <- numeric()
    for(feature_index in 1:length(feature_rank)){
      if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q){
        break
      }else{
        null_feature <- c(null_feature, feature_rank[feature_index])
      }
    }
    MDS_selected_index <- setdiff(feature_rank, null_feature)
  }else{
    MDS_selected_index = NULL
  }

  MDS_selected_index <- MDS_selected_index[order(MDS_selected_index)]

  return(colnames(MS)[MDS_selected_index])
}


#' @title High Dimensional Mediation Analysis using Multiple Data Splitting (MDS)
#' @description The main function of CoxMDS which combines the screening step, data splitting step and estimation step.
#' @param X a vector of the exposure.
#' @param Y a vector of the outcome.
#' @param M a matrix of the mediators.
#' @param COV optional, a matrix of the potential covariates.
#' @param penalty the penalty used in the cox regression.
#' @param intercept binary, whether an intercept should be included in the regression M~X, the default value is TRUE.
#' @param q1 the pre-defined FDR level in the screening step.
#' @param q2 the pre-defined FDR level in the data splitting procedure.
#' @param n_split the pre-defined number of splittings in the MDS.
#'
#' @return A list of the CoxMDS results including the IDs of the screening step, the IDs of the selection step and the estimation for the effects.
#' @export


CoxMDS <- function(X, Y, M, COV, penalty=c('MCP', 'lasso', 'SCAD'), intercept=TRUE, q1=0.2,
                   q2=0.1, n_split = 25){

  n <- nrow(M); p <- ncol(M)

  ##### STEP 1: Preliminary Screening #####

  screen.results <- screen_FDR(X, M, COV, intercept = intercept, fdr=q1)

  MS <- M[, screen.results$ID_screen]

  ##### STEP 2: Selection with Data Splitting #####

  ID_DS <- datasplit_surv(X, Y, MS, COV, penalty = penalty,q=q2, n_split = n_split)

  ID_DS <- screen.results$ID_screen[ID_DS]

  ##### STEP 3: Estimation of Effects #####

  alpha_est <- screen.results$coef[ID_DS]

  DATA <- data.frame(Y = Y, M[, ID_DS], X = X, COV = COV)
  cox_model <- survival::coxph(Y ~ ., data = DATA)

  beta_est <- summary(cox_model)$coefficients[1: length(ID_DS)]

  result <- data.frame(alpha=alpha_est, beta=beta_est, `alpha*beta` = alpha_est*beta_est)

  return(list(first.step=names(screen.results$ID_screen),
              second.step=names(ID_DS),
              estimates=result))
}




########## Useful Functions for MDS ##########

analys <- function(W, q){
  W_abs <- abs(W)
  fdp <- rep(NA, length(W))
  for (i in 1:length(W)) {
    fdp[i] <- sum(W<=-W_abs[i]) / max(1, sum(W>=W_abs[i]))
  }
  id <- which(fdp<=q)
  tau <- min(W_abs[id])
  out <- which(W>=tau)
  return(out)
}









