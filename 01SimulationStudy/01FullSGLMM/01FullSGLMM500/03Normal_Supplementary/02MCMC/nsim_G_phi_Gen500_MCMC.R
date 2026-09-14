################################################################################
################################################################################
# 10192024
# Running Multiple Simulations for Gaussian Data
# MCMC
# Test args
################################################################################
################################################################################

rm(list=ls())
library(cli);library(viridis);library(fields);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm)
library(invgamma);library(pROC);library(Matrix);library(emulator)
library(argparser);library(tidyverse);


#Distributing our work to different CPU
args <- arg_parser('jin_test') %>%
  add_argument('--task', type = 'integer', help = 'Task number.') %>%
  parse_args()

phi_test<-c(1,2,3,4) #it corresponds to c(0.1, 0.3, 0.5, 0.7)
nsim_test<-1:50 #should change this to 50 later
settings<-expand.grid(phi_test=phi_test,nsim_test=nsim_test)
phiselect<-settings$phi_test[args$task]
sim<-settings$nsim_test[args$task]


file_name <- paste0("../01DataGeneration/G_phi", phiselect, "_Gen500_", sim, ".RData")
load(file_name)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e5, 1)  # Randomly select a base seed
randomseed <- base_seed + sim + phiselect * 1000
set.seed(randomseed)



################################################################################
################################################################################
## MCMC 
################################################################################
################################################################################
## Gibbs sampling 
n.iteration=100000
## Prior Distribution
## Prior Distribution For Beta
prior_mu_beta<-c(1,1); prior_sigma_beta<-diag(2) 
cholB<-chol(prior_sigma_beta)
inv_sigma_beta<-chol2inv(cholB)
## Prior Distribution For sigma2
alpha1=2; beta1=2
## Distribution of Posterior Distribution
## Making vacant matrix for MCMC samples for (Beta1, Beta2, Sigma2,phi)
parameterMatrix<-matrix(NA,nrow=n.iteration, ncol=4)
parameterMatrix[1,]<-c(0,0,1,0.5)
z_minus_xb<-as.numeric(Z_trainingdata-X_trainingdata%*%parameterMatrix[1,c(1,2)])
## Prior for the phi
propSD<-0.1
lb<-0
ub<-sqrt(1)
cholR<-chol(R_phi_trainingdata)
invR_phi_trainingdata<-chol2inv(cholR)
detR<- -(sum(log(diag(cholR))))
pt<-proc.time()
for ( i in 2: n.iteration){
  if(i%%1000==0) {print(i)}
  #Updating Beta
  chol_covmat_beta<-chol((1/parameterMatrix[i-1,3])*t(X_trainingdata)%*%invR_phi_trainingdata%*%(X_trainingdata)+inv_sigma_beta)
  covmat_beta<-chol2inv(chol_covmat_beta)
  mu_beta<-covmat_beta%*%((1/parameterMatrix[i-1,3])*t(X_trainingdata)%*%invR_phi_trainingdata%*%(Z_trainingdata)+inv_sigma_beta%*%prior_mu_beta)
  parameterMatrix[i,c(1,2)]<-mu_beta+t(chol(covmat_beta))%*%rnorm(2)
  z_minus_xb<-as.numeric(Z_trainingdata-X_trainingdata%*%parameterMatrix[i,c(1,2)])
  
  #Updating Sigma2
  alpha1_tilde=alpha1+(0.8*n)/2
  beta1_tilde=beta1+as.numeric(t(z_minus_xb)%*%invR_phi_trainingdata%*%(z_minus_xb))*1/2
  parameterMatrix[i,3]<-1/rgamma(n=1,shape=alpha1_tilde,rate=beta1_tilde)
  
  #update Phi
  proposed_phi<-rnorm(1, mean=parameterMatrix[i-1,4], sd=propSD)
  if
  (proposed_phi<lb | proposed_phi>ub){
    alpha<-(-1e10) }
  else {
    R_prop_phi<-exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/proposed_phi)
    chol_R_prop_phi<-chol(R_prop_phi)
    invR_prop_phi<-chol2inv(chol_R_prop_phi)
    det_R_prop_phi<- -(sum(log(diag(chol_R_prop_phi))))
    
    num<- (-0.8*n/2)*log(parameterMatrix[i,3])+det_R_prop_phi-(1/(2*parameterMatrix[i,3]))*(t(Z_trainingdata-X_trainingdata%*%parameterMatrix[i,c(1,2)])%*%invR_prop_phi%*%(Z_trainingdata-X_trainingdata%*%parameterMatrix[i,c(1,2)]))
    
    den<- (-0.8*n/2)*log(parameterMatrix[i,3])+detR-(1/(2*parameterMatrix[i,3]))*(t(Z_trainingdata-X_trainingdata%*%parameterMatrix[i,c(1,2)])%*%invR_phi_trainingdata%*%(Z_trainingdata-X_trainingdata%*%parameterMatrix[i,c(1,2)]))
    
    alpha<-min(0,num-den)}
  uniff<-runif(1)
  
  if(alpha>log(uniff)){
    invR_phi_trainingdata<-invR_prop_phi
    detR<-det_R_prop_phi
    parameterMatrix[i,4]<-proposed_phi
  }
  else {
    parameterMatrix[i,4]<-parameterMatrix[i-1,4]
  }
}
ptFinal<-proc.time()-pt
ptFinal

################################################################################
################################################################################
## Checking Acceptance rate 
## acceptance rate for Beta1
acceptRate_beta1<-1-mean(parameterMatrix[-n.iteration,1]==parameterMatrix[-1,1])
## acceptance rate for Beta2
acceptRate_beta2<-1-mean(parameterMatrix[-n.iteration,2]==parameterMatrix[-1,2])
## acceptance rate for sigma2
acceptRate_sigma2<-1-mean(parameterMatrix[-n.iteration,3]==parameterMatrix[-1,3])
## Acceptance Rate for Phi
acceptRate_phi<-1-mean(parameterMatrix[-n.iteration,4]==parameterMatrix[-1,4])
## burninperiod
burninperiod<-ceiling(n.iteration*0.1)


################################################################################
################################################################################
## Prediction_MCMC
################################################################################
################################################################################
MCMCsigma2samples <- parameterMatrix[, 3]
MCMCphisamples <- parameterMatrix[, 4]


## Prediction_MCMC
Z_prediction_MCMC<-matrix(rep(NA,(0.2*n)*n.iteration), nrow=(0.2*n), ncol=n.iteration)

for(i in 1:n.iteration){
  if (i %% 1000 == 0) { print(i) }
  Sigma_11 <- MCMCsigma2samples[i] * exp(-distMat[c(1:(0.8 * n)), c(1:(0.8 * n))] / MCMCphisamples[i])
  chol_inv_Sigma_11 <- chol(Sigma_11)
  inv_Sigma_11 <- chol2inv(chol_inv_Sigma_11)
  Sigma_21 <- MCMCsigma2samples[i] * exp(-distMat[c((0.8 * n + 1):n), c(1:(0.8 * n))] / MCMCphisamples[i])
  Z_prediction_MCMC[,i]<-X_testdata%*%parameterMatrix[i,c(1,2)]+Sigma_21%*%inv_Sigma_11%*%(Z_trainingdata-X_trainingdata%*%parameterMatrix[i,c(1,2)])
}
Z_prediction_MCMC
Z_prediction_MCMC_mean<-apply(Z_prediction_MCMC,1,mean)
Z_RMSPE_MCMC<-sqrt(mean((Z_testdata-Z_prediction_MCMC_mean)^2))


################################################################################
################################################################################
save(parameterMatrix, ptFinal, Z_RMSPE_MCMC,Z_prediction_MCMC,Z_prediction_MCMC_mean,
     acceptRate_beta1, acceptRate_beta2, acceptRate_sigma2, acceptRate_phi,n.iteration,burninperiod,
     n.iteration,
     file = paste0("G_phi",phiselect,"_Gen500_",sim,"_MCMC.RData"))
