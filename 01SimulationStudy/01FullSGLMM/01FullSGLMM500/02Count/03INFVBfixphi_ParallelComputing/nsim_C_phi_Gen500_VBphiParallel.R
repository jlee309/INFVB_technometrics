################################################################################
################################################################################
# 10142024
# Running Multiple Simulations for Count Data
# VB_Parallel : Fixing Phi
################################################################################
################################################################################
rm(list=ls())
library(cli);library(viridis);library(fields);library(pgdraw);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm);
library(pROC);library(Matrix);library(emulator);
library(snow);library(doParallel);library(foreach);library(parallel);

detectCores() # Checks the number of cores on computer
nprocs<-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # You have the flexibility to select the number of cores for parallel computing in this line. In our analysis, we utilized 30 cores, as specified in the slurm file.
this_cluster <- parallel::makeCluster(nprocs, type="PSOCK") # Initiates the Cluster
doParallel::registerDoParallel(this_cluster) # Register the cluster

## INFVB 
# thetaC={Gamma(Beta,W,Sigma2)} : Estimating Gamma, Sigma2 and Phi
# thetaD={phi} : Discretizng this Phi; Do not Estimate phi. 
################################################################################
################################################################################
## Preliminaries
library(invgamma) 
library(mvtnorm)


n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed



for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name <- paste0("../01DataGeneration/C_phi", phiselect, "_Gen500_", sim, ".RData")
    load(file_name)
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)


p<-length(beta)
iter<-1000 # (Note) Different from n.iteration in previous MCMC

## VB Beta Prior (Same as MCMC prior)
VB_mu_beta_prior<-c(0,0); VB_sigma_beta_prior<-diag(2) 
chol_VB_sigma_beta_prior<-chol(VB_sigma_beta_prior)
VB_inv_Sigma_beta_prior<-chol2inv(chol_VB_sigma_beta_prior)

## VB Sigma2 Prior
VB_iter<-matrix(NA, nrow=iter, ncol=2) 
VB_iter[1,]<-c(1,1) #Initial values
VB_alpha_sigma<-2 ; VB_beta_sigma<-2
VB_E_sigma2Inv<-VB_iter[1,1]/VB_iter[1,2]


## Disciticizing ThetaD:{phi} 
iter.thetaphi=1000 #Number of Discriticiaing Theta 
thetaphi<-seq(0.001,sqrt(2),length.out=iter.thetaphi)


## ELBOvector 
ELBOvector<-ELBOjthcal<-numeric()


## Final Saving Matrix for Gamma(Beta, W) and sigma2
thetaC_Meangamma<-matrix(NA, nrow=iter.thetaphi, ncol=((0.8*n)+p))
dim(thetaC_Meangamma)
thetaC_covMatgamma<-matrix(NA, nrow=iter.thetaphi, ncol=(0.8*n+p))
dim(thetaC_covMatgamma)
thetaC_sigma2<-matrix(NA, nrow=iter.thetaphi, ncol=2)

## VB_invR_phi for training data
R_phi_trainingdata<-exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/phi)
VB_invR_phi_trainingdata<-chol2inv(chol(R_phi_trainingdata))

## Making X_tilde(X+Identity)
X_tilde_trainingdata<-cbind(X_trainingdata,diag(0.8*n))
dim(X_tilde_trainingdata)

## Making Gamma(beta+W) Matrix
VB_MeanbetaMat<-matrix(NA, nrow=iter, ncol=2) 
VB_MeanWMat<-matrix(NA, nrow=iter, ncol=(0.8*n)) #Only for training W
VB_MeangammaS<-cbind(VB_MeanbetaMat, VB_MeanWMat)
VB_MeangammaS[1,]<-rep(1,(0.8*n+p))
dim(VB_MeangammaS)
VB_covMatgammaS<-matrix(NA, nrow=iter, ncol=(0.8*n+p))
dim(VB_covMatgammaS)
VB_covMatgammaS[1,]<-rep(1,(0.8*n+p))

#Making Prior inverse covariance matrix of Gamma. 
zeroMat2n<-matrix(0, nrow=p, ncol=(0.8*n)) 
zeroMatn2<-matrix(0, nrow=(0.8*n), ncol=p)
invgammaMatmaking1<-cbind(VB_inv_Sigma_beta_prior,zeroMat2n)
invgammaMatmaking2<-cbind(zeroMatn2,(1/sigma2)*VB_invR_phi_trainingdata)
invgammaMat<-rbind(invgammaMatmaking1,invgammaMatmaking2)
dim(invgammaMat)


################################################################################
################################################################################
VB_pt<-proc.time()
outputMat<-foreach::foreach(j=1:iter.thetaphi,.combine = "cbind",
                            .packages = c("mvtnorm","invgamma","fields")) %dopar% {
                              
                              # Updating phi
                              R_phi_trainingdata<-exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/thetaphi[j]) 
                              cholR_trainingdata<-chol(R_phi_trainingdata)
                              VB_invR_phi_trainingdata<-chol2inv(cholR_trainingdata) 
                              ELBOjthcal<-numeric()
                              ELBOjthcal[1]<-10
                              
                              k=1
                              bar<-77 #Arbitrary number 
                              
                              while(bar>rep(1e-3)) {
                                
                                k<- k+1
                                print(k)
                                
                                # Updatng Discreticized Phi and E[1/(sigma^2)] for invgammaMat and gammaMat
                                invgammaMatmaking2<-cbind(zeroMatn2,VB_E_sigma2Inv*VB_invR_phi_trainingdata)
                                invgammaMat<-rbind(invgammaMatmaking1,invgammaMatmaking2)
                                cholinvgammaMat<-chol(invgammaMat)
                                gammaMat<-chol2inv(cholinvgammaMat)
                                cholgammaMat<-chol(gammaMat)
                                
                                # gamma (Using Optim function) by using Laplace approximation 
                                f_optimgamma<-function(gamma){t(obs_trainingdata)%*%(X_tilde_trainingdata%*%gamma)-t(rep(1,(0.8*n)))%*%exp(X_tilde_trainingdata%*%gamma)-1/2*t(gamma)%*%invgammaMat%*%(gamma)-sum(log(diag(cholgammaMat))) }
                                
                                grG<-function(gamma){t(X_tilde_trainingdata)%*%obs_trainingdata-t(X_tilde_trainingdata)%*%exp(X_tilde_trainingdata%*%gamma)-invgammaMat%*%(gamma)}
                                
                                gamma_afteroptim<-optim(par=c(rep(0,(0.8*n+p))),method="BFGS",
                                                        control=list(fnscale=-1),
                                                        fn=f_optimgamma,gr=grG, hessian = TRUE)
                                
                                VB_Meangamma<-gamma_afteroptim$par
                                VB_covMatgamma<-(solve(-gamma_afteroptim$hessian))
                                VB_MeangammaS[k,]<-as.numeric(VB_Meangamma)
                                VB_covMatgammaS[k,]<-diag(VB_covMatgamma)
                                VB_SigmaW<-VB_covMatgamma[c((p+1):(0.8*n+p)),c((p+1):(0.8*n+p))] 
                                
                                #Sigma2 (Using conjugacy between sigma2 and W)
                                VB_iter[k,1]<-VB_alpha_sigma+(0.8*n)/2
                                VB_iter[k,2]<-VB_beta_sigma+0.5*(t(VB_MeangammaS[k,c((p+1):(0.8*n+p))])%*%VB_invR_phi_trainingdata%*%(VB_MeangammaS[k,c((p+1):(0.8*n+p))])+sum(diag(VB_invR_phi_trainingdata%*%VB_SigmaW)))
                                VB_E_sigma2Inv<-VB_iter[k,1]/VB_iter[k,2]
                                
                                #ELBO preparation
                                VB_covMatgamma_chol<-chol(VB_covMatgamma)
                                
                                # ELBO 
                                ELBO1<- t(obs_trainingdata)%*%(X_tilde_trainingdata%*%VB_Meangamma)-t(rep(1,(0.8*n)))%*%exp(X_tilde_trainingdata%*%VB_Meangamma+(1/2)*diag(X_tilde_trainingdata%*%VB_covMatgamma%*%t(X_tilde_trainingdata)))
                                
                                ELBO2<- -sum(log(diag(cholgammaMat)))-(1/2)*t(VB_Meangamma)%*%(invgammaMat)%*%(VB_Meangamma)-(1/2)*sum(diag(invgammaMat%*%VB_covMatgamma))
                                
                                ELBO3<- (VB_alpha_sigma-VB_iter[k,1])*(digamma(VB_iter[k,1])-log(VB_iter[k,2]))+(VB_iter[k,2]-VB_beta_sigma)*VB_E_sigma2Inv
                                
                                ELBO4<- sum(log(diag(VB_covMatgamma_chol)))-VB_iter[k,1]*log(VB_iter[k,2])
                                
                                
                                ELBOjthcal[k] <- ELBO1+ELBO2+ELBO3+ELBO4
                                
                                #Making ELBO for the Stopping Criteria
                                bar<-abs(ELBOjthcal[k]-ELBOjthcal[k-1])
                              }
                              
                              return(c(ELBOjthcal[k],VB_MeangammaS[k,],VB_iter[k,],VB_covMatgammaS[k,]))
                            }
VB_ptFinal<-proc.time()-VB_pt
VB_ptFinal


################################################################################
################################################################################
ELBOvector<-outputMat[1,]
# Exponential ELBO
weights<-exp(ELBOvector-max(ELBOvector))
ELBOvectorweights<-weights/sum(weights)

################################################################################
################################################################################
thetaC_Meangamma<-matrix(NA, nrow=iter.thetaphi, ncol=((0.8*n)+p))
dim(thetaC_Meangamma)
thetaC_sigma2<-matrix(NA, nrow=iter.thetaphi, ncol=2)
dim(thetaC_sigma2)
thetaC_covMatgamma<-matrix(NA, nrow=iter.thetaphi, ncol=((0.8*n)+p))
dim(thetaC_covMatgamma)


for ( i in 1:iter.thetaphi){
  thetaC_Meangamma[i,]<-outputMat[2:((0.8*n)+p+1),i]
}

for ( i in 1:iter.thetaphi){
  thetaC_sigma2[i,]<-outputMat[((0.8*n)+p+2):((0.8*n)+p+3),i]
}

for ( i in 1:iter.thetaphi){
  thetaC_covMatgamma[i,]<-outputMat[((0.8*n)+p+4):(1+((0.8*n)+p)+2+((0.8*n)+p)),i]
}



################################################################################
################################################################################
# Making new Density according to ELBOvector rates
################################################################################
################################################################################
# Beta1 Weighted Density 
par(mfrow=c(1,1))
matDensityBeta1<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqb1<-seq(-2,2, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesBeta1 <-dnorm(x=xSeqb1, mean=thetaC_Meangamma[k,1], sd=sqrt(thetaC_covMatgamma[k,1]))
  matDensityBeta1[,k]<-densitiesBeta1
}
dim(matDensityBeta1)
newdensityBeta1<-matDensityBeta1%*%ELBOvectorweights
#plot(x=xSeqb1, y=newdensityBeta1, typ="l", pch=16, col="red",main="Beta1")

# Beta2 Weighted Density 
matDensityBeta2<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqb2<-seq(-2,2, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesBeta2 <-dnorm(x=xSeqb2, mean=thetaC_Meangamma[k,p], sd=sqrt(thetaC_covMatgamma[k,2]))
  matDensityBeta2[,k]<-densitiesBeta2
}
dim(matDensityBeta2)
newdensityBeta2<-matDensityBeta2%*%ELBOvectorweights
dim(newdensityBeta2)
#plot(x=xSeqb2, y=newdensityBeta2, typ="l", pch=16, col="red", xlim=c(0,2), main="Beta2")


#Sigma2 Weighted Density
library(invgamma)
matDensitySigma2<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqsigma2<-seq(0.1,6, length.out=iter.thetaphi)
matDensitySigma2<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesSigma2 <-dinvgamma(x=xSeqsigma2, shape=thetaC_sigma2[k,1], rate=thetaC_sigma2[k,2])
  matDensitySigma2[,k]<-densitiesSigma2
}
dim(matDensitySigma2)
newdensitySigma2<-matDensitySigma2%*%ELBOvectorweights
dim(newdensitySigma2)
#plot(x=xSeqsigma2, y=newdensitySigma2, typ="l", pch=16, col="red", main="Sigma2")


################################################################################
################################################################################
## W(Spatial random Effects)
################################################################################
################################################################################
# W3 Weighted Density (Wqth) 
q=3
matDensityW3<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqW3<-seq(-3,3, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesW3 <-dnorm(x=xSeqW3, mean=thetaC_Meangamma[k,(p+q)], sd=sqrt(thetaC_covMatgamma[k,p+q]))
  matDensityW3[,k]<-densitiesW3
}
dim(matDensityW3)
newdensityW3<-matDensityW3%*%ELBOvectorweights

# W10 Weighted Density (Wqth) 
q=10
matDensityW10<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqW10<-seq(-3,3, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesW10 <-dnorm(x=xSeqW10, mean=thetaC_Meangamma[k,(p+q)], sd=sqrt(thetaC_covMatgamma[k,p+q]))
  matDensityW10[,k]<-densitiesW10
}
dim(matDensityW10)
newdensityW10<-matDensityW10%*%ELBOvectorweights

# W25 Weighted Density (Wqth) 
q=25
matDensityW25<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqW25<-seq(-3,3, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesW25 <-dnorm(x=xSeqW25, mean=thetaC_Meangamma[k,(p+q)], sd=sqrt(thetaC_covMatgamma[k,p+q]))
  matDensityW25[,k]<-densitiesW25
}
dim(matDensityW25)
newdensityW25<-matDensityW25%*%ELBOvectorweights

# W40 Weighted Density (Wqth) 
q=40
matDensityW40<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqW40<-seq(-3,3, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesW40 <-dnorm(x=xSeqW40, mean=thetaC_Meangamma[k,(p+q)], sd=sqrt(thetaC_covMatgamma[k,p+q]))
  matDensityW40[,k]<-densitiesW40
}
dim(matDensityW40)
newdensityW40<-matDensityW40%*%ELBOvectorweights




################################################################################
################################################################################
## Prediction_VB
################################################################################
################################################################################
lambda_prediction_VB<-matrix(rep(NA,(0.2*n)*iter.thetaphi), nrow=(0.2*n), ncol=iter.thetaphi)
eta_trainingdata_VB<-matrix(rep(NA,(0.8*n)*iter.thetaphi), nrow=(0.8*n), ncol=iter.thetaphi)
eta_testdata_VB<-matrix(rep(NA,(0.2*n)*iter.thetaphi), nrow=(0.2*n), ncol=iter.thetaphi)


for(i in 1:iter.thetaphi){
  if(i%%100==0) {print(i)}
  Sigma_11<-(thetaC_sigma2[i,2]/(thetaC_sigma2[i,1]-1))*exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/thetaphi[i])
  inv_Sigma_11<-solve(Sigma_11)
  Sigma_21<-(thetaC_sigma2[i,2]/(thetaC_sigma2[i,1]-1))*exp(-distMat[c((0.8*n+1):(n)),c(1:(0.8*n))]/thetaphi[i])
  Sigma_22<-(thetaC_sigma2[i,2]/(thetaC_sigma2[i,1]-1))*exp(-distMat[c((0.8*n+1):(n)),c((0.8*n+1):(n))]/thetaphi[i])
  Sigma_12<-(thetaC_sigma2[i,2]/(thetaC_sigma2[i,1]-1))*exp(-distMat[c(1:(0.8*n)),c((0.8*n+1):(n))]/thetaphi[i])
  
  eta_trainingdata_VB[,i]<-X_trainingdata%*%thetaC_Meangamma[i,c(1,2)]+thetaC_Meangamma[i,c((p+1):(0.8*n+p))]
  eta_testdata_VB[,i]<-X_testdata%*%thetaC_Meangamma[i,c(1,2)]+Sigma_21%*%(inv_Sigma_11%*%thetaC_Meangamma[i,c((p+1):(0.8*n+p))])
  lambda_prediction_VB[,i]<-exp(eta_testdata_VB[,i])
}

lambda_prediction_VB<-lambda_prediction_VB%*%ELBOvectorweights
lambda_RMSPE_VB<-sqrt(mean((lambda_testdata-lambda_prediction_VB)^2))
obs_RMSPE_VB<-sqrt(mean((obs_testdata-lambda_prediction_VB)^2))


################################################################################
################################################################################
# save(obs_RMSPE_VB,VB_ptFinal,ELBOvectorweights,newdensityBeta1, newdensityBeta2,
#      newdensitySigma2,thetaphi,iter.thetaphi,lambda_prediction_VB,lambda_RMSPE_VB,
#      newdensityW3, newdensityW10, newdensityW25, newdensityW40,file="C_phi01_Gen500_Results_VBphiParallel.RData")


save(obs_RMSPE_VB,VB_ptFinal,newdensityBeta1, newdensityBeta2,
     newdensitySigma2,thetaphi,iter.thetaphi,ELBOvectorweights,
     newdensityW3, newdensityW10, newdensityW25, newdensityW40,
     lambda_prediction_VB,lambda_RMSPE_VB, 
     thetaC_covMatgamma,thetaC_Meangamma,thetaC_sigma2,
     file=paste0("C_phi",phiselect,"_Gen500_",sim,"_VBphiParallel.RData")) 

  }
}

