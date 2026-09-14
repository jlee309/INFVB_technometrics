################################################################################
################################################################################
# MFVB for Count Data using 100 Eigen Basis functions 
# Rcpp Functions are not used

################################################################################
################################################################################
rm(list=ls())
library(fields) ; library(MASS) ; library(Matrix)
library(invgamma)  ; library(mvtnorm) ; library(emulator)



################################################################################
################################################################################
# Number of simulations
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name1 <- paste0("../../01DataGeneration/phi",phiselect,"_Gen25k_",sim,"EigenBasis.RData")
    file_name2 <-paste0("../../01DataGeneration/phi",phiselect,"_Gen25k_",sim,"SpatialData.RData")
    load(file_name1)
    load(file_name2)
    
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    
################################################################################
################################################################################
# Fast Self Made Functions
quad.diag.New<-function (M, x) { # Diagonal elements of a quadratic form diag(t(x)%*%M%*%x)
  colSums(crossprod(M, x) * x)
}
quad.form.new<-function(M,x){ # t(x)%*%M%*%x
  crossprod(crossprod(M, x), x)
}

computeTrace<-function(A,B){  #sum(diag(A%*B))
  sum(colSums(t(A)*B))
}

# Set Basis Functions
################################################################################
################################################################################
numofbasis<-100 # Choose the first 100 eigen basis functions
M<-refEigen_trainingdata[,1:numofbasis]
M_CV<-refEigen_testdata[,1:numofbasis]
dim(M);dim(M_CV)
################################################################################
################################################################################
################################################################################
## MFVB(Laplace Approximations)
################################################################################
################################################################################
################################################################################
################################################################################
## Preliminaries
## Initial Settings for MFVB
p<-length(beta)
iter<-1000

## Beta
VB_mean_beta<-c(0,0)
VB_covMat_beta<-100*diag(2)
VB_inv_Sigma_beta_prior<-solve(VB_covMat_beta)

## Sigma2
## Sigma saving matrix 
MFVB_iter<-matrix(NA, nrow=iter, ncol=2) 
MFVB_iter[1,]<-c(1,1) #Initial values

## Prior Distribution For Sigma2 (Inverse Gamma) 
VB_alpha_sigma<-0.1 ; VB_beta_sigma<-0.1
VB_E_sigma2Inv<-MFVB_iter[1,1]/MFVB_iter[1,2]


## Prior Distribution For Delta 
VB_mu_delta<-c(rep(0,numofbasis)) ;VB_sigma_delta<-VB_E_sigma2Inv*diag(numofbasis)
VB_cholD<-chol(VB_sigma_delta)
VB_inv_sigma_delta<-chol2inv(VB_cholD)


#Making X_tilde(X_trainingdata+refRadial_trainingdata)
X_tilde<-cbind(X_trainingdata,M)
tX_tilde<-t(X_tilde)
dim(X_tilde)

#Making Gamma(beta+delta) Matrix
VB_MeanbetaMat<-matrix(NA, nrow=iter, ncol=p) 
dim(VB_MeanbetaMat)
VB_mean_deltaMat<-matrix(NA,nrow=iter,ncol=(numofbasis)) 
dim(VB_mean_deltaMat)

#First Line
VB_MeanbetaMat[1,]<-c(rep(1,p))
VB_mean_deltaMat[1,]<-c(rep(0,numofbasis))

#Making MFVB_Meangamma Matrix
MFVB_Meangamma<-cbind(VB_MeanbetaMat, VB_mean_deltaMat)
dim(MFVB_Meangamma)

#Making covMatgamma Matrix 
MFVB_covMatgamma<-matrix(NA, nrow=iter, ncol=numofbasis+p)
dim(MFVB_covMatgamma)
MFVB_covMatgamma[1,]<-rep(1,numofbasis+p)
zeroMatpn<-matrix(0, nrow=p, ncol=numofbasis) 
zeroMatnp<-matrix(0, nrow=numofbasis, ncol=p)


#Prior
invgammaMatmaking1<-cbind(VB_inv_Sigma_beta_prior,zeroMatpn)
dim(invgammaMatmaking1)
invgammaMatmaking2<-cbind(zeroMatnp,VB_inv_sigma_delta)
dim(invgammaMatmaking2)
invgammaMat<-rbind(invgammaMatmaking1,invgammaMatmaking2)
dim(invgammaMat)


#Expensive Operations
tXObs<-as.numeric(t(X_tilde)%*%obsPois_trainingdata)

#ELBO vector
ELBOjthcal<-numeric()
ELBOjthcal[1]<-10


basisInd<-p+(1:numofbasis)
betaInd<-1:p

MFVB_pt<-proc.time()
for (k in 2:iter){
  print(k)
  # Updatng invgammaMat 
  VB_inv_sigma_delta<-diag(numofbasis) 
  invgammaMatmaking2<-cbind(zeroMatnp,Diagonal(x=rep(VB_E_sigma2Inv,numofbasis)))
  invgammaMat<-rbind(invgammaMatmaking1,invgammaMatmaking2)
  
  # gamma (Using Optim function)
  f_optimgamma<-function(gamma){as.numeric(tXObs%*%gamma-sum(exp(X_tilde%*%gamma))-1/2*t(gamma)%*%invgammaMat%*%(gamma)) }
  grG<-function(gamma){as.numeric(tXObs-t(X_tilde)%*%exp(X_tilde%*%gamma)-invgammaMat%*%(gamma))}
  gamma_afteroptim<-optim(par=c(rep(0,(numofbasis+2))),method="BFGS",
                          control=list(fnscale=-1),
                          fn=f_optimgamma,gr=grG, hessian = TRUE)
  VB_mean_gamma<-gamma_afteroptim$par
  VB_covMat_gamma<-(solve(-gamma_afteroptim$hessian))
  MFVB_Meangamma[k,]<-as.numeric(VB_mean_gamma)
  MFVB_covMatgamma[k,]<-diag(VB_covMat_gamma)
  Sigmadelta<-VB_covMat_gamma[basisInd,basisInd]
  
  #Sigma2
  MFVB_iter[k,1]<-VB_alpha_sigma+numofbasis/2
  MFVB_iter[k,2]<-VB_beta_sigma+0.5*(t(MFVB_Meangamma[k,basisInd])%*%(MFVB_Meangamma[k,basisInd])+sum(diag(Sigmadelta)))
  VB_E_sigma2Inv<-MFVB_iter[k,1]/MFVB_iter[k,2]
  
################################################################################
  # Stopping Criteria (By using difference between ELBO)
  # ELBO preparation
  MFVB_covMatgamma_chol<-chol(VB_covMat_gamma)

  ELBO1<- tXObs%*%VB_mean_gamma-t(rep(1,length(obsPois_trainingdata)))%*%exp(X_tilde%*%VB_mean_gamma+(1/2)*quad.diag.New(M=VB_covMat_gamma,x = tX_tilde)) 
  ELBO2<- -(1/2)*t(VB_mean_gamma)%*%(invgammaMat%*%VB_mean_gamma)-(1/2)*computeTrace(invgammaMat,VB_covMat_gamma)
  ELBO3<- -VB_beta_sigma*(MFVB_iter[k,1]/MFVB_iter[k,2])+sum(log(diag(MFVB_covMatgamma_chol)))-MFVB_iter[k,1]*log(MFVB_iter[k,2])

  
  ELBOjthcal[k] <- ELBO1+ELBO2+ELBO3
  ELBObar<-abs(ELBOjthcal[k]-ELBOjthcal[k-1])< 1e-4
  if(all(ELBObar)){MFVB_Meangamma<-MFVB_Meangamma[1:k,];MFVB_covMatgamma<-MFVB_covMatgamma[1:k,];MFVB_iter<-MFVB_iter[1:k,];break} 

}
MFVB_ptFinal<-proc.time()-MFVB_pt
MFVB_ptFinal


################################################################################
################################################################################
## Prediction obs_testdata 
################################################################################
################################################################################
#Prediction_VB
j<-nrow(MFVB_Meangamma)
lambda_prediction_MFVB<-exp(X_testdata%*%MFVB_Meangamma[j,c(1,2)]
                          +M_CV%*%MFVB_Meangamma[j,basisInd])
lambda_RMSPE_MFVB<-sqrt(mean((lambda_prediction_MFVB-lambda_testdata)^2))

exp_RMSPE_VB<-vector("numeric")
for(j in 1:nrow(MFVB_Meangamma)){
  exp_prediction_VB<-exp(X_testdata%*%MFVB_Meangamma[j,c(1,2)]+M_CV%*%MFVB_Meangamma[j,basisInd])
  exp_RMSPE_VB[j]<-sqrt(mean((exp_prediction_VB-obsPois_testdata)^2))
}


obs_RMSPE_MFVB<-sqrt(mean((lambda_prediction_MFVB-obsPois_testdata)^2))

save(lambda_prediction_MFVB , lambda_RMSPE_MFVB , obs_RMSPE_MFVB ,MFVB_ptFinal , 
     MFVB_Meangamma , MFVB_covMatgamma ,MFVB_iter , 
     file=paste0("100Basis_C_phi",phiselect,"_Gen25k_",sim,"_HMFVB.RData")) 
  }
}
