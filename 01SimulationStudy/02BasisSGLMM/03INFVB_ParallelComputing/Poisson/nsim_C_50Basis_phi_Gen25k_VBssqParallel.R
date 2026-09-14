################################################################################
################################################################################
# INFVB for Count Data using 50 Basis functions 
# when sigma2 is discretized (Parallelization Computing)  
# Rcpp functions are not used 
################################################################################
################################################################################
rm(list=ls())
getwd()
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
    


# Set Basis Functions
################################################################################
################################################################################
numofbasis<-50 
M<-refEigen_trainingdata[,1:numofbasis]
M_CV<-refEigen_testdata[,1:numofbasis]
dim(M);dim(M_CV)
################################################################################
################################################################################
################################################################################
## VB(Laplace Approximations)
################################################################################
################################################################################
################################################################################
################################################################################
## Preliminaries
## Initial Settings for VB
p<-length(beta)
iter<-1000

## Beta
VB_mean_beta<-c(0,0)
VB_covMat_beta<-100*diag(2)
VB_inv_Sigma_beta_prior<-solve(VB_covMat_beta)


## Prior Distribution For Sigma2 (Inverse Gamma) 
VB_alpha_sigma<-0.1 ; VB_beta_sigma<-0.1
# VB_E_sigma2Inv<-VB_iter[1,1]/VB_iter[1,2]


## Prior Distribution For Delta 
VB_mu_delta<-c(rep(0,numofbasis)) ;VB_sigma_delta<-1*diag(numofbasis)
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

#Making VB_Meangamma Matrix
VB_Meangamma<-cbind(VB_MeanbetaMat, VB_mean_deltaMat)
dim(VB_Meangamma)

#Making covMatgamma Matrix 
VB_covMatgamma<-matrix(NA, nrow=iter, ncol=numofbasis+p)
dim(VB_covMatgamma)
VB_covMatgamma[1,]<-rep(1,numofbasis+p)
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


# ELBOvector 
ELBOvector<-ELBOjthcal<-numeric()

basisInd<-p+(1:numofbasis)
betaInd<-1:p

################################################################################
# INFVB Setting 
# Discretizing ThetaD:{sigma2} 
iter.thetaphi=1000 #Number of Discretizing ThetaD 
x_sigma2=seq(0.001,2000,length.out=iter.thetaphi)
theta_sigma2=dinvgamma(x=x_sigma2,shape=VB_alpha_sigma, rate=VB_beta_sigma)
#plot(x_sigma2,theta_sigma2,col="blue")

# Final Saving Matrix for V(Beta, W) and sigma2
thetaC_Meangamma<-matrix(NA, nrow=iter.thetaphi, ncol=(numofbasis+p))
dim(thetaC_Meangamma)
thetaC_covMatgamma<-matrix(NA, nrow=iter.thetaphi, ncol=(numofbasis+p))
dim(thetaC_covMatgamma)

################################################################################
################################################################################
VB_pt<-proc.time()
outputMat<-foreach::foreach(j=1:iter.thetaphi,.combine = "cbind",.packages =c("mvtnorm","invgamma","fields","Matrix")) %dopar% {
  
  # Updatng invgammaMat 
  VB_inv_sigma_delta<-diag(numofbasis) 
  invgammaMatmaking2<-cbind(zeroMatnp,Diagonal(x=rep((1/x_sigma2[j]),numofbasis)))
  invgammaMat<-rbind(invgammaMatmaking1,invgammaMatmaking2)
  cholinvgammaMat<-chol(invgammaMat)
  
  ELBO3<-log(theta_sigma2[j])
  ELBOjthcal[1]<-10
  k=1
  bar<-77 #Arbitrary number 
  
  while(bar>exp(1e-3)) {
    k<-k+1
    print(k)
  # gamma (Using Optim function)
  f_optimgamma<-function(gamma){as.numeric(tXObs%*%gamma-sum(exp(X_tilde%*%gamma))-1/2*t(gamma)%*%invgammaMat%*%(gamma)) }
  grG<-function(gamma){as.numeric(tXObs-t(X_tilde)%*%exp(X_tilde%*%gamma)-invgammaMat%*%(gamma))}
  gamma_afteroptim<-optim(par=c(rep(0,(numofbasis+2))),method="BFGS",
                          control=list(fnscale=-1),
                          fn=f_optimgamma,gr=grG, hessian = TRUE)
  VB_mean_gamma<-gamma_afteroptim$par
  VB_covMat_gamma<-(solve(-gamma_afteroptim$hessian))
  VB_Meangamma[k,]<-as.numeric(VB_mean_gamma)
  VB_covMatgamma[k,]<-diag(VB_covMat_gamma)

################################################################################
  # Stopping Criteria (By using difference between ELBO)
  # ELBO calculation 
  VB_covMatgamma_chol<-chol(VB_covMat_gamma)
  
  ELBO1<- tXObs%*%VB_mean_gamma-t(rep(1,length(obsPois_trainingdata)))%*%exp(X_tilde%*%VB_mean_gamma+(1/2)*quad.diag.New(M=VB_covMat_gamma,x = tX_tilde)) 
  ELBO2<- -(1/2)*t(VB_mean_gamma)%*%(invgammaMat%*%VB_mean_gamma)-(1/2)*computeTrace(invgammaMat,VB_covMat_gamma)
  ELBO4<- sum(log(diag(VB_covMatgamma_chol)))+sum(log(diag(cholinvgammaMat)))
  
  ELBOjthcal[k] <- ELBO1+ELBO2+ELBO3+ELBO4
  
  #Making ELBO for the Stopping Criteria
  bar<-abs(ELBOjthcal[k]-ELBOjthcal[k-1])
  }
  
  ELBOvector<-c(ELBOvector,ELBOjthcal[k])
  thetaC_Meangamma[j,]<-VB_Meangamma[k,]
  thetaC_covMatgamma[j,]<-VB_covMatgamma[k,]
  
  return(c(ELBOjthcal[k],VB_Meangamma[k,],VB_covMatgamma[k,]))
}
VB_ptFinal<-proc.time()-VB_pt
VB_ptFinal

dim(outputMat)
ELBOvector<-outputMat[1,]
# Exponential ELBO
weights<-exp(ELBOvector-max(ELBOvector))
ELBOvectorweights<-weights/sum(weights)



################################################################################
################################################################################
for ( i in 1:iter.thetaphi){
  thetaC_Meangamma[i,]<-outputMat[2:((numofbasis)+p+1),i]
}
for ( i in 1:iter.thetaphi){
  thetaC_covMatgamma[i,]<-outputMat[((numofbasis)+p+2):(1+p+numofbasis+p+numofbasis),i]
}


########################################################################
########################################################################
# Making new Density according to ELBOvector rates
########################################################################
########################################################################
# Beta1 Weighted Density 
par(mfrow=c(1,1))
matDensityBeta1<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqb1<-seq(-4,4, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesBeta1 <-dnorm(x=xSeqb1, mean=thetaC_Meangamma[k,1], sd=sqrt(thetaC_covMatgamma[k,1]))
  matDensityBeta1[,k]<-densitiesBeta1
}
dim(matDensityBeta1)
newdensityBeta1<-matDensityBeta1%*%ELBOvectorweights
#plot(x=xSeqb1, y=newdensityBeta1, typ="l", pch=16, col="red",main="Beta1")

# Beta2 Weighted Density 
matDensityBeta2<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqb2<-seq(-4,4, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesBeta2 <-dnorm(x=xSeqb2, mean=thetaC_Meangamma[k,p], sd=sqrt(thetaC_covMatgamma[k,2]))
  matDensityBeta2[,k]<-densitiesBeta2
}
dim(matDensityBeta2)
newdensityBeta2<-matDensityBeta2%*%ELBOvectorweights
dim(newdensityBeta2)
#plot(x=xSeqb2, y=newdensityBeta2, typ="l", pch=16, col="red", xlim=c(0,2), main="Beta2")


############################################################
## delta
############################################################
# delta3 Weighted Density
q=3
matDensitydelta3<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqdelta3<-seq(-10,10, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesdelta3 <-dnorm(x=xSeqdelta3, mean=thetaC_Meangamma[k,(p+q)], sd=sqrt(thetaC_covMatgamma[k,p+q]))
  matDensitydelta3[,k]<-densitiesdelta3
}
dim(matDensitydelta3)
newdensitydelta3<-matDensitydelta3%*%ELBOvectorweights


################################################################################
################################################################################
## Prediction obs_testdata 
################################################################################
################################################################################
#Prediction_VB
#j<-nrow(VB_Meangamma)
lambda_prediction_VB<-matrix(rep(NA,(0.2*n)*iter.thetaphi), nrow=(0.2*n), ncol=iter.thetaphi)
eta_testdata_VB<-matrix(rep(NA,(0.2*n)*iter.thetaphi), nrow=(0.2*n), ncol=iter.thetaphi)
eta_trainingdata_VB<-matrix(rep(NA,(0.8*n)*iter.thetaphi), nrow=(0.8*n), ncol=iter.thetaphi)

dim(lambda_prediction_VB)
for(i in 1:iter.thetaphi){
  if(i%%100==0) {print(i)}
  eta_trainingdata_VB[,i]<-X_trainingdata%*%thetaC_Meangamma[i,c(1,2)]+M%*%thetaC_Meangamma[i,c((1+p):(p+numofbasis))]
  eta_testdata_VB[,i]<-X_testdata%*%thetaC_Meangamma[i,c(1,2)]+M_CV%*%thetaC_Meangamma[i,c((1+p):(p+numofbasis))]
  lambda_prediction_VB[,i]<-exp(eta_testdata_VB[,i])
}
lambda_prediction_VB<-lambda_prediction_VB%*%ELBOvectorweights
lambda_RMSPE_VB<-sqrt(mean((lambda_prediction_VB-lambda_testdata)^2))
lambda_RMSPE_VB


obs_RMSPE_VB<-sqrt(mean((lambda_prediction_VB-obsPois_testdata)^2))


save(lambda_prediction_VB , lambda_RMSPE_VB , obs_RMSPE_VB ,VB_ptFinal , ELBOvectorweights,
     xSeqb1, xSeqb2, xSeqdelta3,newdensityBeta1,newdensityBeta2,newdensitydelta3,x_sigma2,iter.thetaphi,
     thetaC_Meangamma , thetaC_covMatgamma, 
     eta_testdata_VB,eta_trainingdata_VB,numofbasis,p,
file=paste0("50Basis_C_phi",phiselect,"_Gen25k_",sim,"_VBssqParallel.RData")) 
  }
}
