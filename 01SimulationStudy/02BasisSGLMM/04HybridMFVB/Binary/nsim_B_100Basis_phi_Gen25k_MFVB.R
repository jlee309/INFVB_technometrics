################################################################################
################################################################################
# HMFVB for Binary Data using 100 Eigen Basis functions 
# Rcpp functions are used 
# 10262024
# Using 
################################################################################
################################################################################
rm(list=ls())
library(fields) ; library(MASS) ; library(Matrix)
library(invgamma)  ; library(mvtnorm) ; library(emulator)
library(Rcpp) # Load package 'Rcpp'
library(RcppEigen) # Load package 'RcppArmadillo'

# sourceCpp("matsourceEigen.cpp")
sourceCpp("matsourceEigen.cpp", cacheDir = "/scratch/jlee309/rcpp-cache/Basis_INFVB_phi")



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
numofbasis<-100 # Choose the first 100 Eigen basis functions
M<-refEigen_trainingdata[,1:numofbasis]
M_CV<-refEigen_testdata[,1:numofbasis]
################################################################################
################################################################################
# Variational Bayes (MFVB)
## Preliminaries
## Initial Settings for MFVB
p<-length(beta)
n<-length(obsBin_trainingdata)
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


## Making C matrix 
dim(refEigen_trainingdata)
C<-cbind(X_trainingdata,M)
tC<-t(C)
dim(C)


#Making V(beta+W) Matrix
VB_MeanbetaMat<-matrix(NA, nrow=iter, ncol=p) 
dim(VB_MeanbetaMat)
VB_mean_deltaMat<-matrix(NA,nrow=iter,ncol=(numofbasis)) 
dim(VB_mean_deltaMat)


#First Line
VB_MeanbetaMat[1,]<-c(rep(1,p))
VB_mean_deltaMat[1,]<-c(rep(0,numofbasis))


#Making MFVB_MeanV Matrix
MFVB_MeanV<-cbind(VB_MeanbetaMat, VB_mean_deltaMat)
dim(MFVB_MeanV)


#Making covMatV Matrix 
MFVB_covMatV<-matrix(NA, nrow=iter, ncol=((numofbasis)+p))
dim(MFVB_covMatV)
MFVB_covMatV[1,]<-rep(1,((numofbasis)+p))
zeroMatpn<-matrix(0, nrow=p, ncol=(numofbasis)) 
zeroMatnp<-matrix(0, nrow=(numofbasis), ncol=p)


#Prior
invVMatmaking1<-cbind(VB_inv_Sigma_beta_prior,zeroMatpn)
dim(invVMatmaking1)
invVMatmaking2<-cbind(zeroMatnp,VB_E_sigma2Inv*VB_inv_sigma_delta)
dim(invVMatmaking2)
invVMat<-rbind(invVMatmaking1,invVMatmaking2)
dim(invVMat)

#ELBO vector
ELBOjthcal<-numeric()
ELBOjthcal[1]<-10

# Observations
obsMinHalf<-obsBin_trainingdata-0.5
tCobsMinHalf<-as.numeric(tC%*%obsMinHalf)

################################################################################
################################################################################
#Making temporary lambda to get xi (Using Rcpp Functions)
#Method A-1(Without Using Rcpp)
#system.time({xi<-sqrt(quad.diag.New(M=(MFVB_covMatV[1,]+MFVB_MeanV[1,]%*%t(MFVB_MeanV[1,])),x = tC))})

#Method A-2(With Using Rcpp) #Method A-2 is much faster than Method A-1
system.time({
  foo<-(MFVB_covMatV[1,]+MFVB_MeanV[1,]%*%t(MFVB_MeanV[1,]))
  xi<-(sqrt(myfun3(foo,tC))) 
})

lambda<- -tanh(xi/2)/(4*xi)
lambda<-Diagonal(x = lambda)

#Method B-1(Without Using Rcpp) 
#system.time({tCLambdaC_1<-quad.form.new(M=lambda,x=C)}) 

#Method B-2(With Using Rcpp) #Method B-2 is much faster than B-1
system.time({
  baz<-as.matrix(lambda%*%C)
  tCLambdaC<-myfun2(C,baz)}) 


################################################################################
################################################################################
#MFVB Preparation

VB_covMatchol<-chol(invVMat -2*tCLambdaC)
VB_covMatVing<-chol2inv(VB_covMatchol)
VB_covMatVing<-as.matrix(VB_covMatVing)
MFVB_covMatV[1,]<- diag(VB_covMatVing)
MFVB_MeanV[1,] <- as.numeric(VB_covMatVing%*%tCobsMinHalf)


m<-numofbasis
basisInd<-(1:m)+p
betaInd<-1:p
################################################################################
################################################################################
#MFVB Start
MFVB_pt<-proc.time()
for (k in 2:iter){
  print(k)
  
  #Sigma2
  Sigmaw<-VB_covMatVing[basisInd,basisInd]
  MFVB_iter[k,1]<-VB_alpha_sigma+(numofbasis)/2
  MFVB_iter[k,2]<-VB_beta_sigma+0.5*(t(MFVB_MeanV[k-1,basisInd])%*%(MFVB_MeanV[k-1,basisInd])+sum(diag(Sigmaw)))
  VB_E_sigma2Inv<-MFVB_iter[k,1]/MFVB_iter[k,2]
  
  #Update xi- Computational Bottlenecks with Matrix -Matrix quadratic forms
  foo<-(VB_covMatVing+MFVB_MeanV[k-1,]%*%t(MFVB_MeanV[k-1,]))
  xi<-sqrt(myfun3(foo,tC)) 
  lambda<- -tanh(xi/2)/(4*xi)
  sum_psi_xi<-sum(xi*0.5-log(1+exp(xi))+xi*(tanh(xi/2)/(4)))
  lambda<-Diagonal(x = lambda)
  baz<-as.matrix(lambda%*%C)
  tCLambdaC<-myfun2(C,baz)
  
  # Updatng invVMat 
  invVMatmaking2<-cbind(zeroMatnp,Diagonal(x=rep(VB_E_sigma2Inv,m)))
  invVMat<-rbind(invVMatmaking1,invVMatmaking2)
  VB_covMatchol<-chol(invVMat -2*tCLambdaC)
  VB_covMatVing<-chol2inv(VB_covMatchol)
  VB_covMatVing<-as.matrix(VB_covMatVing)
  MFVB_covMatV[k,]<- diag(VB_covMatVing)
  MFVB_MeanV[k,] <- as.numeric(VB_covMatVing%*%tCobsMinHalf)
  
  #############################################################
  # # Stopping Criteria (By using difference between ELBO)
  #ELBO preparation
  ELBO1<-sum(quad.form.new(M =tCLambdaC , x= MFVB_MeanV[k,]),
             sum(diag(tCLambdaC%*%VB_covMatVing)),
             as.numeric(-0.5*quad.form.new(M = invVMat, x= MFVB_MeanV[k,])),
             -0.5*sum(diag(invVMat%*%VB_covMatVing)), 
             sum_psi_xi)

  ELBO2<- tCobsMinHalf%*%MFVB_MeanV[k,]
  ELBO3<- -VB_beta_sigma*(MFVB_iter[k,1]/MFVB_iter[k,2])-MFVB_iter[k,1]*log(MFVB_iter[k,2])
  ELBO4<- -sum(log(diag(VB_covMatchol))) #same as 1/2*log(det(VB_covMatVing))
  ELBOjthcal[k] <- ELBO1+ELBO2+ELBO3+ELBO4
  
  ELBObar<-abs(ELBOjthcal[k]-ELBOjthcal[k-1])<rep(1e-4)
  #print(ELBOjthcal)
  
  if(all(ELBObar)){MFVB_MeanV<-MFVB_MeanV[1:k,];MFVB_covMatV<-MFVB_covMatV[1:k,];MFVB_iter<-MFVB_iter[1:k,];break}
}
MFVB_ptFinal<-proc.time()-MFVB_pt
MFVB_ptFinal
################################################################################
################################################################################
## Cross Validation
################################################################################
################################################################################
j<-nrow(MFVB_MeanV)
#Prediction_VB
eta_trainingdata_MFVB<-X_trainingdata%*%MFVB_MeanV[j,c(1,2)]+M%*%MFVB_MeanV[j,basisInd]
eta_testdata_MFVB<-X_testdata%*%MFVB_MeanV[j,c(1,2)]+M_CV%*%MFVB_MeanV[j,basisInd]
exp_prediction_MFVB<-exp(eta_testdata_MFVB)/(1+exp(eta_testdata_MFVB))
exp_RMSPE_MFVB<-sqrt(mean((exp_prediction_MFVB-obsBin_testdata)^2))

library(pROC)
rocBinary <- roc(obsBin_testdata, as.numeric(exp_prediction_MFVB))
aucVal_MFVB<-auc(rocBinary)

save(exp_prediction_MFVB , exp_RMSPE_MFVB , aucVal_MFVB ,MFVB_ptFinal , 
     MFVB_MeanV , MFVB_covMatV ,MFVB_iter , 
     eta_testdata_MFVB, eta_trainingdata_MFVB, 
    file=paste0("100Basis_B_phi",phiselect,"_Gen25k_",sim,"_HMFVB.RData")) 
  }
}
