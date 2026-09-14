################################################################################
################################################################################
# 11072024 
# Basis SGLMM Stan file for Binary 
################################################################################
################################################################################
rm(list=ls())
library(fields) ; library(mvtnorm)
library(rstan) ; library(nimble);library(pROC)
source("sharedFunctions.R")
source("batchmeans.r")


# Number of simulations
n_simulations <- 1 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed

for(phiselect in 2:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name1 <- paste0("../../01DataGeneration/phi",phiselect,"_Gen25k_",sim,"EigenBasis.RData")
    file_name2 <-paste0("../../01DataGeneration/phi",phiselect,"_Gen25k_",sim,"SpatialData.RData")
    load(file_name1)
    load(file_name2)
    
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)


# Set Basis Functions
################################################################################
################################################################################
numofbasis<-20 
M<-refEigen_trainingdata[,1:numofbasis]
M_CV<-refEigen_testdata[,1:numofbasis]
################################################################################
################################################################################
# Variational Bayes (VB)
## Preliminaries
## Initial Settings for VB
k<-ncol(X_trainingdata)
p<-length(beta)
nn<-length(obsBin_trainingdata)

################################################################################
# Stan Setup
################################################################################
# Fit Model Using Stan
iter=1000
gamma_samples<-stan(file="BasisSGLMM.stan", 
                    data = list(N=nn,
                                K=k,
                                P=numofbasis,
                                y=obsBin_trainingdata,
                                X=X_trainingdata,
                                M=M),
                    pars = c("betas","deltas","ssq"),
                    iter=iter,
                    warmup = floor(iter/10),
                    chains = 1)

#Extract Pertinent Information from output
totTime<-sum(get_elapsed_time(gamma_samples)) # Model Fitting Time
parMat<-cbind(rstan::extract(gamma_samples, c('betas'))[[1]], # Extract samples for Model Parameters
              rstan::extract(gamma_samples, c('ssq'))[[1]])
deltaMat<-rstan::extract(gamma_samples, c('deltas'))[[1]] # Extract samples for basis coefficients

# print(gamma_samples)
# stan_trace(gamma_samples,pars = c('betas','ssq'))
################################################################################################
# Save Data
################################################################################################
# save(parMat,deltaMat,totTime,
#      file="mcmc_Stan_small.RData") # Save Data
################################################################################################
# Summarize Samples
# load("../output/mcmc_Stan_small.RData")
################################################################################################
summaryMat<-list()
summaryMat[[1]]<-round(summaryFunction(parMat,
                                       totTime=totTime),3)
summaryMat[[2]]<-round(summaryFunction(deltaMat,
                                       totTime=totTime),3)
# summaryMat[[1]] # Summary Table
# apply(summaryMat[[2]],1,mean) # Mean results for the basis coefficients

################################################################################################
# Compute Mean Squared Prediction Error and AUC
################################################################################################
#Dr.Lee's way to do this.
mcmcDat<-list()
mcmcDat[[1]]<-parMat; colnames(mcmcDat[[1]])<-c("beta1","beta2","ssq")
mcmcDat[[2]]<-deltaMat
cvSummary<-cvFunction.gamma(mcmcDat=mcmcDat,mBase=M_CV,
                            XMatCV=X_testdata,obsCV=obsBin_testdata)
# cvSummary[1]
# cvSummary[2] 
# print(cvSummary[[1]]) #MSPE
# totTime/60

#Jin's way 
# dim(parMat)
# dim(deltaMat)
# n

#Prediction_MCMC
lengthPredSamples<-iter-iter/10
exp_prediction<-matrix(rep(NA,(0.2*n)*lengthPredSamples), nrow=(0.2*n), ncol=lengthPredSamples)
dim(exp_prediction)
eta_testdata_MCMC_Stan <- matrix(rep(NA, (0.2 * n) * lengthPredSamples), nrow = (0.2 * n), ncol = lengthPredSamples)
eta_trainingdata_MCMC_Stan <- matrix(rep(NA, (0.8 * n) * lengthPredSamples), nrow = (0.8 * n), ncol = lengthPredSamples)

dim(eta_testdata_MCMC_Stan)
dim(eta_trainingdata_MCMC_Stan)

for(i in 1:lengthPredSamples){
  eta_trainingdata_MCMC_Stan[, i]<-X_trainingdata%*%parMat[i,c(1,2)]+M%*%deltaMat[i,] 
  eta_testdata_MCMC_Stan[,i]<-X_testdata%*%parMat[i,c(1,2)]+M_CV%*%deltaMat[i,]
  exp_prediction[,i] <-as.numeric(exp(eta_testdata_MCMC_Stan[,i])/(1+exp(eta_testdata_MCMC_Stan[,i])))
  
}
exp_prediction<-apply(exp_prediction,1,mean) 
exp_RMSPE_MCMC_Stan<-sqrt(mean((exp_prediction-obsBin_testdata)^2))
# exp_RMSPE_MCMC_Stan


rocBinary <- roc(obsBin_testdata, exp_prediction)
# Get the full AUC
aucVal_Stan<-auc(rocBinary)
aucVal_Stan

################################################################################################
# Save Data
################################################################################################
save(summaryMat,parMat,deltaMat,totTime,cvSummary,
     aucVal_Stan, exp_RMSPE_MCMC_Stan, 
     eta_testdata_MCMC_Stan, eta_trainingdata_MCMC_Stan,
    file = paste0("Basis_B_phi",phiselect,"_Gen25k_",sim,"_MCMC_Stan.RData"))
  }
}
