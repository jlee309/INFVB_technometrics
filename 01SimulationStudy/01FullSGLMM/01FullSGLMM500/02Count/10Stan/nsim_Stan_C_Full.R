################################################################################
################################################################################
# 11092024 
# Full SGLMM Stan file for Count Data
################################################################################
################################################################################
rm(list=ls())
library(fields) ; library(mvtnorm)
library(rstan) ; library(nimble);library(pROC)
source("sharedFunctions.R")
source("batchmeans.r")
n_simulations <- 1 
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  
base_seed <- sample.int(1e6, 1) 
for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name <- paste0("../01DataGeneration/C_phi", phiselect, "_Gen500_", sim, ".RData")
    load(file_name)
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
k<-ncol(X_trainingdata)
p<-length(beta)
nn<-length(obs_trainingdata)
nnt<-length(obs_testdata)
M<-diag(nn)
M_CV<-diag(nnt)
iter=1000
distMat_M<-distMat[c(1:(0.8*n)),c(1:(0.8*n))]
dim(distMat_M)
gamma_samples<-stan(file="FullSGLMM_C.stan", 
                    data = list(N=nn,
                                K=k,
                                P=nn,
                                y=obs_trainingdata,
                                X=X_trainingdata,
                                M=M,
                                distMat=distMat_M),
                    pars = c("betas","omegas","ssq","phi"),
                    iter=iter,
                    warmup = floor(iter/10),
                    chains = 1)
totTime<-sum(get_elapsed_time(gamma_samples))
parMat<-cbind(rstan::extract(gamma_samples, c('betas'))[[1]], 
              rstan::extract(gamma_samples, c('phi'))[[1]],
              rstan::extract(gamma_samples, c('ssq'))[[1]])
omegasMat<-rstan::extract(gamma_samples, c('omegas'))[[1]]
stan_trace(gamma_samples,pars = c('betas','ssq','phi'))
summaryMat<-list()
summaryMat[[1]]<-round(summaryFunction(parMat,
                                       totTime=totTime),3)
summaryMat[[2]]<-round(summaryFunction(omegasMat,
                                       totTime=totTime),3)
mcmcDat<-list()
mcmcDat[[1]]<-parMat; colnames(mcmcDat[[1]])<-c("beta1","beta2","phi","ssq")
mcmcDat[[2]]<-omegasMat
lengthPredSamples<-iter-iter/10
MCMCsigma2samples <- parMat[, 4]
MCMCphisamples <- parMat[, 3]
lambda_prediction<-matrix(rep(NA,(0.2*n)*lengthPredSamples), nrow=(0.2*n), ncol=lengthPredSamples)
dim(lambda_prediction)
eta_testdata_MCMC_Stan <- matrix(rep(NA, (0.2 * n) * lengthPredSamples), nrow = (0.2 * n), ncol = lengthPredSamples)
eta_trainingdata_MCMC_Stan <- matrix(rep(NA, (0.8 * n) * lengthPredSamples), nrow = (0.8 * n), ncol = lengthPredSamples)
CovMat_prediction_MCMC_Stan<-matrix(rep(NA,(0.2*n)^2*lengthPredSamples), nrow=(0.2*n)^2, ncol=lengthPredSamples)
dim(eta_testdata_MCMC_Stan)
dim(eta_trainingdata_MCMC_Stan)
for(i in 1:lengthPredSamples){
  if (i %% 10== 0) { print(i) }
  Sigma_11 <- MCMCsigma2samples[i] * exp(-distMat[c(1:(0.8 * n)), c(1:(0.8 * n))] / MCMCphisamples[i])
  chol_inv_Sigma_11 <- chol(Sigma_11)
  inv_Sigma_11 <- chol2inv(chol_inv_Sigma_11)
  Sigma_21 <- MCMCsigma2samples[i] * exp(-distMat[c((0.8 * n + 1):n), c(1:(0.8 * n))] / MCMCphisamples[i])
  Sigma_22<-MCMCphisamples[i]*exp(-distMat[c((0.8*n+1):(n)),c((0.8*n+1):(n))]/MCMCphisamples[i])
  Sigma_12<-MCMCphisamples[i]*exp(-distMat[c(1:(0.8*n)),c((0.8*n+1):(n))]/MCMCphisamples[i])
  eta_trainingdata_MCMC_Stan[, i]<-X_trainingdata%*%parMat[i,c(1,2)]+omegasMat[i,] 
  eta_testdata_MCMC_Stan[,i]<-X_testdata%*%parMat[i,c(1,2)]+Sigma_21 %*% (inv_Sigma_11 %*%omegasMat[i,])
  lambda_prediction[,i] <-as.numeric(exp(eta_testdata_MCMC_Stan[,i]))
  # Prediction 
  CovMat_prediction_MCMC_Stan[,i]<-as.numeric(Sigma_22-Sigma_21%*%inv_Sigma_11%*%Sigma_12)
}
lambda_prediction<-apply(lambda_prediction,1,mean) 
lambda_RMSPE_Stan<-sqrt(mean((lambda_prediction-lambda_testdata)^2))
lambda_RMSPE_Stan
obs_RMSPE_MCMC_Stan<-sqrt(mean((lambda_prediction-obs_testdata)^2))
obs_RMSPE_MCMC_Stan
save(summaryMat,parMat,omegasMat,totTime,
     obs_RMSPE_MCMC_Stan, lambda_RMSPE_Stan,
     eta_testdata_MCMC_Stan, eta_trainingdata_MCMC_Stan,CovMat_prediction_MCMC_Stan,
    file = paste0("C_phi",phiselect,"_Gen500_",sim,"_MCMC_Stan.RData"))}}
