################################################################################
################################################################################
# INFVB for Binary Data using 50 Eigen Basis functions
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


file_name1 <- paste0("../../01DataGeneration/phi",phiselect,"_Gen25k_",sim,"EigenBasis.RData")
file_name2 <-paste0("../../01DataGeneration/phi",phiselect,"_Gen25k_",sim,"SpatialData.RData")
load(file_name1)
load(file_name2)


set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed
randomseed <- base_seed + sim + phiselect * 1000
set.seed(randomseed)


library(fields) ; library(MASS) ; library(Matrix)
################################################################################
################################################################################
# Set Basis Functions
################################################################################
################################################################################
numofbasis<-50
M<-refEigen_trainingdata[,1:numofbasis]
M_CV<-refEigen_testdata[,1:numofbasis]
dim(M);dim(M_CV)
################################################################################
################################################################################
## MCMC(Metropolis Hastings within GIBBS sampling)
################################################################################
################################################################################
n.iteration=100000
## Prior Distribution
## Prior Distribution For Beta
mu_beta<-c(0,0); sigma_beta<-100*diag(2) 
cholB<-chol(sigma_beta)
inv_sigma_beta<-chol2inv(cholB)
## Prior Distribution For Delta 
mu_delta<-c(rep(0,numofbasis)) ; sigma_delta<-sigma2*diag(numofbasis) 
cholD<-chol(sigma_delta)
inv_sigma_delta<-chol2inv(cholD)
dim(inv_sigma_delta)
## Prior Distribution For sigma2
alpha1=0.1; beta1=0.1 


## Distribution of Posterior Distribution
## Making vacant matrix for MCMC samples for (Beta1, Beta2, Sigma2)
parameterMatrix<-matrix(NA,nrow=n.iteration, ncol=3)
parameterMatrix[1,]<-c(0,0,1)
dim(parameterMatrix)
## Making vacant matrix for Spatial_Basis-delta coefficient
deltaMatrix<-matrix(NA,nrow=n.iteration,ncol=numofbasis)
deltaMatrix[1,]<-rnorm(numofbasis)
dim(deltaMatrix)
library(MASS)
pt<-proc.time()
for ( i in 2:n.iteration){
  
  if(i%%100==0) {print(i)}
  
  #Estimating Betas(MH within Gibbs)
  proposed_beta<-rnorm(n=2, mean=parameterMatrix[i-1,c(1,2)], sd= 0.06324555) #Proposal Distribution
  probProp<-as.numeric(1/(1+exp(-X_trainingdata%*%proposed_beta-M%*%deltaMatrix[i-1,])))
  num<-(sum(dbinom(x=obsBin_trainingdata,size=1,prob=probProp,log = TRUE)))-1/2*t(proposed_beta)%*%(inv_sigma_beta%*%proposed_beta)
  probPrev<-as.numeric(1/(1+exp(-X_trainingdata%*%parameterMatrix[i-1,c(1,2)]-M%*%deltaMatrix[i-1,])))
  den<-(sum(dbinom(x=obsBin_trainingdata,size=1,prob=probPrev, log = TRUE)))-1/2*t(parameterMatrix[i-1,c(1,2)])%*%(inv_sigma_beta%*%parameterMatrix[i-1,c(1,2)])
  alpha<-min(0,num-den)
  if(alpha>log(runif(1))){
    parameterMatrix[i,c(1,2)]<-proposed_beta
  }else {
    parameterMatrix[i,c(1,2)]<-parameterMatrix[i-1,c(1,2)]
  }
  
  #Updating Sigma2(By Gibbs Sampling)
  alpha1_tilde=alpha1+numofbasis/2
  beta1_tilde=beta1+1/2*t(deltaMatrix[i-1,])%*%(inv_sigma_delta%*%(deltaMatrix[i-1,]))
  parameterMatrix[i,3]<-1/rgamma(n=1,shape=alpha1_tilde ,rate=beta1_tilde)
  
  #Estimating delta(MH within Gibbs),(Every single component) 
  ##### I made changes here - All at once
  proposed_delta<-rnorm(n=numofbasis, mean=deltaMatrix[i-1,], sd=1.1)
  probProp<-as.numeric(1/(1+exp(-X_trainingdata%*%parameterMatrix[i,c(1,2)]-M%*%proposed_delta)))
  num<-(sum(dbinom(x=obsBin_trainingdata,size=1,prob=probProp,log = TRUE)))-1/(2*parameterMatrix[i,3])*t(proposed_delta)%*%(inv_sigma_delta%*%(proposed_delta))
  probPrev<-as.numeric(1/(1+exp(-X_trainingdata%*%parameterMatrix[i,c(1,2)]-M%*%deltaMatrix[i-1,])))
  den<-(sum(dbinom(x=obsBin_trainingdata,size=1,prob=probPrev,log = TRUE)))-1/(2*parameterMatrix[i,3])*t(deltaMatrix[i-1,])%*%(inv_sigma_delta%*%(deltaMatrix[i-1,]))
  
  alpha<-min(0,num-den)
  if(alpha>log(runif(1))){
    deltaMatrix[i,]<-proposed_delta
  }else{
    deltaMatrix[i,]<-deltaMatrix[i-1,]
  }
}
ptFinal<-proc.time()-pt
ptFinal

################################################################################
################################################################################
##  Acceptance rates for Beta1, Beta2, Sigma2, Delta
################################################################################
################################################################################
## acceptance rate for Beta1
acceptRate_beta1<-1-mean(parameterMatrix[-n.iteration,1]==parameterMatrix[-1,1])
acceptRate_beta1
## acceptance rate for Beta2
acceptRate_beta2<-1-mean(parameterMatrix[-n.iteration,2]==parameterMatrix[-1,2])
acceptRate_beta2
## acceptance rate for Sigma2 (#Gibbs sampling/ Acceptance rate should be 1)
acceptRate_marginal_sigma2_delta<-1-mean(parameterMatrix[-n.iteration,3]==parameterMatrix[-1,3])
acceptRate_marginal_sigma2_delta
## Checking Delta Matrix #Arbitrary pick q=5
q=5
acceptRate_marginal_delta<-1-mean(deltaMatrix[-n.iteration,q]==deltaMatrix[-1,q])
acceptRate_marginal_delta 
################################################################################
################################################################################
################################################################################
## Prediction obs_testdata 
################################################################################
################################################################################
#Prediction_MCMC
lengthPredSamples<-10000
predSequence<-floor(seq(ceiling(n.iteration*0.2), n.iteration, length.out=lengthPredSamples))
exp_prediction<-matrix(rep(NA,(0.2*n)*lengthPredSamples), nrow=(0.2*n), ncol=lengthPredSamples)
dim(exp_prediction)
eta_testdata_MCMC <- matrix(rep(NA, (0.2 * n) * lengthPredSamples), nrow = (0.2 * n), ncol = lengthPredSamples)
eta_trainingdata_MCMC <- matrix(rep(NA, (0.8 * n) * lengthPredSamples), nrow = (0.8 * n), ncol = lengthPredSamples)

for(i in 1:lengthPredSamples){
  useInd<-predSequence[i]
  
  eta_trainingdata_MCMC[, i]<-X_trainingdata%*%parameterMatrix[useInd,c(1,2)]+M%*%deltaMatrix[useInd,] 
  eta_testdata_MCMC[,i]<-X_testdata%*%parameterMatrix[useInd,c(1,2)]+M_CV%*%deltaMatrix[useInd,]
  exp_prediction[,i] <-as.numeric(exp(eta_testdata_MCMC[,i])/(1+exp(eta_testdata_MCMC[,i])))

  }
exp_prediction<-apply(exp_prediction,1,mean) 
exp_RMSPE_MCMC<-sqrt(mean((exp_prediction-obsBin_testdata)^2))
exp_RMSPE_MCMC



library(pROC)
rocBinary <- roc(obsBin_testdata, exp_prediction)
# Get the full AUC
aucVal<-auc(rocBinary)
aucVal


save(parameterMatrix , deltaMatrix , ptFinal , 
     exp_prediction, aucVal , exp_RMSPE_MCMC , 
     acceptRate_marginal_delta , acceptRate_marginal_sigma2_delta , acceptRate_beta2,  
     acceptRate_beta1,n.iteration, 
     eta_trainingdata_MCMC,eta_testdata_MCMC,
     file = paste0("50Basis_B_phi",phiselect,"_Gen25k_",sim,"_MCMC.RData"))


