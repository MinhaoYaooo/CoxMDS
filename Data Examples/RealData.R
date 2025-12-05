library(survival)
library(corpcor)
library(CoxMKF)
library(CoxMDS)

##### DATA PREPROCESSING #####

## Read the methylation data
meth <- read.table('TCGA.LUNG.sampleMap_HumanMethylation450', header = TRUE)
meth <- na.omit(meth)
meth_t <- t(meth)

M <- meth_t[2:nrow(meth_t),1:ncol(meth_t)]
M <- apply(M, 1, as.numeric)
M <- t(M)
M <- as.data.frame(M)

## Read the clinical data
colnames(M) <- meth_t[1,]
pheno <- read.csv('clinical_754.csv', header = TRUE)
pheno$sampleID <- gsub('-','\\.',pheno$sampleID)
common_id <- intersect(pheno$sampleID, row.names(M))  # 754 common id, the same as pheno

Y <- Surv(pheno$OS, pheno$Death)

X <- pheno$X.smoking..current.smoker_1.0.

M <- M[common_id,]    # sort to match the id

n <- nrow(M)

p <- ncol(M)

## Process the baseline covariates
COV <- pheno[, c('age_at_initial_pathologic_diagnosis',
                 'pathologic_stage',
                 'radiation_therapy',
                 'gender')]

colnames(COV) <- c('age', 'stage', 'radio', 'gender')

COV$stage_num <- NA

for (i in 1:n) {
  if(COV$stage[i] %in% c('Stage I', 'Stage IA', 'Stage IB')){
    COV$stage_num[i] <- 1
  }else if(COV$stage[i] %in% c('Stage II', 'Stage IIA', 'Stage IIB')){
    COV$stage_num[i] <- 2
  }else if(COV$stage[i] %in% c('Stage III', 'Stage IIIA', 'Stage IIIB')){
    COV$stage_num[i] <- 3
  }else{
    COV$stage_num[i] <- 4
  }
}

COV$stage <- COV$stage_num
COV$radio <- as.integer(COV$radio=='YES')
COV$gender <- as.integer(COV$gender=='MALE')
COV <- COV[,-5]

##### RUN COXMDS & COXMKF #####

results.MDS <- CoxMDS(X = X, Y = Y, M = M, COV = COV, penalty = "MCP")

results.MKF <- CoxMKF(X = X, Y = Y, M = M, COV = COV, penalty = "MCP)




