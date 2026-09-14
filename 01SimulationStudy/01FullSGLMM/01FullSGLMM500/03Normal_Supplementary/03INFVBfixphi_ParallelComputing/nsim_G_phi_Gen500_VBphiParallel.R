################################################################################
################################################################################
# 10192024
# Running Multiple Simulations for Gaussian Data
# VB_Parallel : Fixing Phi
# Parallel Computing
################################################################################
################################################################################
rm(list=ls())

library(cli);library(viridis);library(fields);library(pgdraw);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm);
library(pROC);library(Matrix);library(emulator);
library(snow);library(doParallel);library(foreach);library(parallel);

detectCores() # Checks the number of cores on computer
nprocs <-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # You have the flexibility to select the number of cores for parallel computing in this line. In our analysis, we utilized 30 cores, as specified in the slurm file.
this_cluster <- parallel::makeCluster(nprocs, type="PSOCK") # Initiates the Cluster
doParallel::registerDoParallel(this_cluster) # Register the cluster
################################################################################
################################################################################
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name <- paste0("../01DataGeneration/G_phi", phiselect, "_Gen500_", sim, ".RData")
    load(file_name)
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    



################################################################################
################################################################################
## INFVB 
# thetaC={Beta,Sigma2)} : Estimating Beta, Sigma2 and Phi
# thetaD={phi} : Discreticizng this Phi; Do not Estimate phi in INFVB. 
################################################################################
################################################################################
## (1) VB Beta Prior
VB_mu_beta_prior<-c(1,1); VB_sigma_beta_prior<-diag(2)
VB_alpha_sigma2<-2 ; VB_beta_sigma2<-2
## Initial Settings for VB (Initial Value for VB) #Don't be confused
VB_mean_beta<-c(0,0)
VB_covMat_beta<-diag(2)
## Preliminaries
iter=1000
VB_iter<-matrix(NA, nrow=iter, ncol=8)
VB_iter[1,]<-c(VB_mean_beta,VB_covMat_beta,2,2)
## Beta update 
chol_VB_sigma_beta_prior<-chol(VB_sigma_beta_prior)
VB_inv_Sigma_beta_prior<-chol2inv(chol_VB_sigma_beta_prior)
## Sigma2
VB_E_sigma2Inv<-VB_iter[1,7]/VB_iter[1,8]

## (2) ThetaD Part
iter.thetaphi=1000
thetaphi<-seq(0.001,1,length.out=iter.thetaphi)



################################################################################
################################################################################
VB_pt<-proc.time()
outputMat<-foreach::foreach(j=1:iter.thetaphi,.combine = "cbind",
                            .packages = c("mvtnorm","invgamma","fields")) %dopar% {
 
  R_phi<-exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/thetaphi[j])
  cholR<-chol(R_phi)
  VB_invR_phi<-chol2inv(cholR)
  ELBOjthcal<-numeric()
  ELBOjthcal[1]<-10

  k = 1
  bar<-33
  
  while(bar>rep(1e-10) ) { 
    k <- k+1
    
    #Sigma2
    VB_iter[k,7]<- VB_alpha_sigma2+(0.8*n)/2
    exp_sigma2 <- t(Z_trainingdata)%*%VB_invR_phi%*%(Z_trainingdata)-2*t(Z_trainingdata)%*%VB_invR_phi%*%(X_trainingdata)%*%(VB_mean_beta)+t(VB_mean_beta)%*%t(X_trainingdata)%*%VB_invR_phi%*%(X_trainingdata)%*%(VB_mean_beta)+sum(diag(t(X_trainingdata)%*%VB_invR_phi%*%(X_trainingdata)%*%VB_covMat_beta))
    VB_iter[k,8]<-VB_beta_sigma2+0.5*exp_sigma2
    VB_E_sigma2Inv<-VB_iter[k,7]/VB_iter[k,8]
    
    # Normal - Beta
    VB_chol_covMat_beta<-chol(VB_E_sigma2Inv*t(X_trainingdata)%*%VB_invR_phi%*%(X_trainingdata)+VB_inv_Sigma_beta_prior)
    VB_covMat_beta<-chol2inv(VB_chol_covMat_beta)
    VB_mean_beta<-as.numeric(VB_covMat_beta%*%(VB_E_sigma2Inv*t(X_trainingdata)%*%VB_invR_phi%*%(Z_trainingdata)+ VB_inv_Sigma_beta_prior%*%VB_mu_beta_prior))
    VB_iter[k,1:2]<-as.numeric(VB_mean_beta)
    VB_iter[k,3:6]<-as.numeric(VB_covMat_beta)
    
    #ELBO(k) Prepartion
    ELBO1 <- -sum(log(diag(cholR)))
    ELBO2 <- -0.5*t(VB_iter[k,c(1:2)]-VB_mu_beta_prior)%*%VB_inv_Sigma_beta_prior%*%(VB_iter[k,c(1:2)]-VB_mu_beta_prior)
    ELBO3 <- -0.5*sum(diag(VB_inv_Sigma_beta_prior%*%matrix(VB_iter[k,c(3:6)], nrow=2, ncol=2,byrow=TRUE)))
    ELBO4<-sum(log(diag(chol(matrix(VB_iter[k,c(3:6)], nrow=2, ncol=2,byrow=TRUE)))))
    ELBO5<- -VB_iter[k,7]*log(VB_iter[k,8])
    ELBOjthcal[k]<- ELBO1+ELBO2+ELBO3+ELBO4+ELBO5
    
    #Making ELBO for the Stopping Criteria
    bar<-abs(ELBOjthcal[k]-ELBOjthcal[k-1])
  }
  return(c(ELBOjthcal[k],VB_iter[k,]))
  # ELBOvector<-c(ELBOvector,ELBOcurrent)
  # thetaC[j,]<-VB_iter[k,]
}
VB_ptFinal<-proc.time()-VB_pt
VB_ptFinal

dim(outputMat)

################################################################################
################################################################################
ELBOvector<-outputMat[1,]
# Exponential ELBO
weights<-exp(ELBOvector-max(ELBOvector))
ELBOvectorweights<-weights/sum(weights)


################################################################################
################################################################################
## (3) ThetaC Part
thetaC_MeanV<-matrix(NA, nrow=iter.thetaphi, ncol=2 )
dim(thetaC_MeanV)
thetaC_sigma2<-matrix(NA, nrow=iter.thetaphi, ncol=2)
dim(thetaC_sigma2)
thetaC_covMatV<-matrix(NA, nrow=iter.thetaphi, ncol=2)
dim(thetaC_covMatV)


for ( i in 1:iter.thetaphi){
  thetaC_MeanV[i,]<-outputMat[2:3,i]
}

for ( i in 1:iter.thetaphi){
  thetaC_sigma2[i,]<-outputMat[8:9,i]
}

for ( i in 1:iter.thetaphi){
  thetaC_covMatV[i,]<-outputMat[c(4,7),i]
}


################################################################################
################################################################################
# Beta1 Weighted Density 
matDensityBeta1<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqb1<-seq(0,2, length.out=iter.thetaphi)
matDensityBeta1<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesBeta1 <-dnorm(x=xSeqb1, mean=thetaC_MeanV[k,1], sd=sqrt(thetaC_covMatV[k,1]))
  matDensityBeta1[,k]<-densitiesBeta1
}
newdensityBeta1<-matDensityBeta1%*%ELBOvectorweights
################################################################################
################################################################################
# Beta2 Weighted Density 
matDensityBeta2<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqb2<-seq(0,2, length.out=iter.thetaphi)
matDensityBeta2<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesBeta2 <-dnorm(x=xSeqb2, mean=thetaC_MeanV[k,2], sd=sqrt(thetaC_covMatV[k,2]))
  matDensityBeta2[,k]<-densitiesBeta2
}
newdensityBeta2<-matDensityBeta2%*%ELBOvectorweights
################################################################################
################################################################################
#Sigma2 Weighted Density 
matDensitySigma2<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqsigma2<-seq(0.1,6, length.out=iter.thetaphi)
matDensitySigma2<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesSigma2 <-dinvgamma(x=xSeqsigma2, shape=thetaC_sigma2[k,1], rate=thetaC_sigma2[k,2])
  matDensitySigma2[,k]<-densitiesSigma2
}
newdensitySigma2<-matDensitySigma2%*%ELBOvectorweights

################################################################################
################################################################################
## Prediction_VB
################################################################################
################################################################################

Sigma_11<-R_phi_trainingdata<-exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/phi)
chol_inv_Sigma_11<-chol(Sigma_11)
inv_Sigma_11<-chol2inv(chol_inv_Sigma_11)
Sigma_21<-exp(-distMat[c((0.8*n+1):(n)),c(1:(0.8*n))]/phi)

## Prediction_VB
Z_prediction_VB<-matrix(rep(NA,(0.2*n)*iter.thetaphi), nrow=(0.2*n), ncol=iter.thetaphi)
# Beta_VB<-cbind(newdensityBeta1,newdensityBeta2)

for(i in 1:iter.thetaphi){
  if(i%%100==0) {print(i)}
  Sigma_11<-(thetaC_sigma2[i,2]/(thetaC_sigma2[i,1]-1))*exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/thetaphi[i])
  inv_Sigma_11<-solve(Sigma_11)
  Sigma_21<-(thetaC_sigma2[i,2]/(thetaC_sigma2[i,1]-1))*exp(-distMat[c((0.8*n+1):(n)),c(1:(0.8*n))]/thetaphi[i])
  # Sigma_22<-(thetaC_sigma2[i,2]/(thetaC_sigma2[i,1]-1))*exp(-distMat[c((0.8*n+1):(n)),c((0.8*n+1):(n))]/thetaphi[i])
  # Sigma_12<-(thetaC_sigma2[i,2]/(thetaC_sigma2[i,1]-1))*exp(-distMat[c(1:(0.8*n)),c((0.8*n+1):(n))]/thetaphi[i])
  Z_prediction_VB[,i]<-X_testdata%*%thetaC_MeanV[i,c(1:2)]+Sigma_21%*%inv_Sigma_11%*%(Z_trainingdata-X_trainingdata%*%thetaC_MeanV[i,c(1:2)])
}
Z_prediction_VB<-Z_prediction_VB%*%ELBOvectorweights
Z_RMSPE_VB<-sqrt(mean((Z_testdata-Z_prediction_VB)^2))


################################################################################
################################################################################
save(Z_RMSPE_VB,VB_ptFinal, Z_prediction_VB,
     ELBOvectorweights,newdensityBeta1, newdensityBeta2,
     newdensitySigma2,thetaphi,iter.thetaphi,xSeqb1,xSeqb2,xSeqsigma2,
     thetaC_MeanV, thetaC_sigma2, thetaC_covMatV,
     file=paste0("G_phi",phiselect,"_Gen500_",sim,"_VBphiParallel.RData")) 
  }
}

