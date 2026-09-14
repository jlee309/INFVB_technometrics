################################################################################
################################################################################
# 10142024
# Running Multiple Simulations for Binary Data
# VB_Parallel : Fixing Phi
################################################################################
################################################################################
rm(list=ls())

library(cli);library(viridis);library(fields);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm)
library(invgamma);library(pROC);library(Matrix);library(emulator)
library(snow);library(doParallel);library(foreach);library(parallel);

detectCores() # Checks the number of cores on computer
nprocs <-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # You have the flexibility to select the number of cores for parallel computing in this line. In our analysis, we utilized 30 cores, as specified in the slurm file.
this_cluster <- parallel::makeCluster(nprocs, type="PSOCK") # Initiates the Cluster
doParallel::registerDoParallel(this_cluster) # Register the cluster


################################################################################
################################################################################
quad.diag.New<-function (M, x) {
  colSums(crossprod(M, x) * x)
}

quad.form.new<-function(M,x){
  crossprod(crossprod(M, x), x)
}
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
    file_name <- paste0("../01DataGeneration/B_phi", phiselect, "_Gen500_", sim, ".RData")
    load(file_name)
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    
## Preliminaries
## Initial Settings for VB
p<-length(beta)
iter<-1000 # (Note) Different from n.iteration in previous MCMC

## Beta
VB_mean_beta<-c(0,0)
VB_covMat_beta<-diag(2)
VB_inv_Sigma_beta_prior<-solve(VB_covMat_beta)


## Sigma2
## Sigma saving matrix 
VB_iter<-matrix(NA, nrow=iter, ncol=2) 
VB_iter[1,]<-c(1,1) #Initial values

## Prior Distribution For Sigma2 (Inverse Gamma) 
VB_alpha_sigma<-0.1 ; VB_beta_sigma<-0.1
VB_E_sigma2Inv<-VB_iter[1,1]/VB_iter[1,2]

## Discretizing ThetaD:{phi} 
iter.thetaphi=1000 #Number of Discretizing Theta 
thetaphi<-seq(0.001,sqrt(2),length.out=iter.thetaphi)

## ELBOvector 
ELBOvector<-ELBOjthcal<-numeric()
#ELBOjthcal[1]<-10 #Pick Arbitrary Number

## W 
R_phi_trainingdata<-exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/phi)
VB_invR_phi_trainingdata<-chol2inv(chol(R_phi_trainingdata))

## Making X_tilde(=previously C ) matrix 
X_tilde<-cbind(X_trainingdata,diag(0.8*n))
tX_tilde<-t(X_tilde)
dim(X_tilde)

#Making V(beta+W) Matrix
VB_MeanbetaMat<-matrix(NA, nrow=iter, ncol=p) 
dim(VB_MeanbetaMat)
VB_mean_wMat<-matrix(NA,nrow=iter,ncol=(0.8*n)) 
dim(VB_mean_wMat)

#First Line
VB_MeanbetaMat[1,]<-c(rep(1,p))
VB_mean_wMat[1,]<-c(rep(0,(0.8*n)))

#Making VB_MeanV Matrix
VB_MeanV<-cbind(VB_MeanbetaMat, VB_mean_wMat)
dim(VB_MeanV)

#Making covMatV Matrix 
VB_covMatV<-matrix(NA, nrow=iter, ncol=((0.8*n)+p))
VB_covMatV[1,]<-c(rep(1,(0.8*n)+p))
zeroMatpn<-matrix(0, nrow=p, ncol=(0.8*n)) 
zeroMatnp<-matrix(0, nrow=(0.8*n), ncol=p)
dim(zeroMatnp)

#Prior
invVMatmaking1<-cbind(VB_inv_Sigma_beta_prior,zeroMatpn)
dim(invVMatmaking1)
invVMatmaking2<-cbind(zeroMatnp,VB_E_sigma2Inv*VB_invR_phi_trainingdata)
dim(invVMatmaking2)
invVMat<-rbind(invVMatmaking1,invVMatmaking2)
dim(invVMat)


# Observations
obsMinHalf<-obs_trainingdata-0.5

#making temporary lambda to get xi 
foo<-(VB_covMatV[1,]+VB_MeanV[1,]%*%t(VB_MeanV[1,]))
xi<-sqrt(quad.diag.New(M=foo,x = tX_tilde))
lambda<- -tanh(xi/2)/(4*xi)
lambda<-diag(x = lambda)
tXtilLambdaXtil<-quad.form.new(M=lambda,x=X_tilde)


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
                                
                                k<-k+1
                                #print(k)
                                
                                # Updatng invVMat 
                                invVMatmaking2<-cbind(zeroMatnp,VB_E_sigma2Inv*VB_invR_phi_trainingdata)
                                invVMat<-rbind(invVMatmaking1,invVMatmaking2)
                                
                                VB_covMatchol<-chol(invVMat -2*tXtilLambdaXtil)
                                VB_covMatVing<-chol2inv(VB_covMatchol)
                                VB_covMatV[k,]<-diag(VB_covMatVing)
                                VB_MeanV[k,] <- VB_covMatVing%*%(tX_tilde%*%(obsMinHalf)) 
                                
                                #Sigma2
                                Sigmaw<-VB_covMatVing[c((p+1):((0.8*n)+p)),c((p+1):((0.8*n)+p))]
                                VB_iter[k,1]<-VB_alpha_sigma+(0.8*n)/2
                                VB_iter[k,2]<-VB_beta_sigma+0.5*(t(VB_MeanV[k,c((p+1):((0.8*n)+p))])%*%VB_invR_phi_trainingdata%*%(VB_MeanV[k,c((p+1):((0.8*n)+p))])+sum(diag(VB_invR_phi_trainingdata%*%Sigmaw)))
                                VB_E_sigma2Inv<-VB_iter[k,1]/(VB_iter[k,2])
                                
                                #Update xi
                                foo<-(VB_covMatVing+VB_MeanV[k,]%*%t(VB_MeanV[k,]))
                                xi<-sqrt(quad.diag.New(M=foo,x = tX_tilde))
                                lambda<- -tanh(xi/2)/(4*xi)
                                lambda<-diag(x = lambda)
                                tXtilLambdaXtil<-quad.form.new(M=lambda,x=X_tilde)
                                
                                #Updata psi(xi)
                                psixi<-xi/2-log(1+exp(xi))+xi*tanh(xi/2)/4
                                psixi<-sum(psixi)
                                #psixi<-diag(x=psixi)
                                
                                #############################################################
                                # # Stopping Criteria (By using difference between ELBO)
                                #ELBO preparation
                                VB_covMatVchol<-chol(VB_covMatVing)
                                invVMat_chol<-chol(invVMat)
                                VMat<-chol2inv(invVMat_chol)
                                VMat_chol<-chol(VMat)
                                
                                ELBO1<-t(VB_MeanV[k,])%*%tX_tilde%*%lambda%*%(X_tilde)%*%VB_MeanV[k,]+sum(diag(tX_tilde%*%(lambda%*%(X_tilde))%*%VB_covMatVing))-0.5*t(VB_MeanV[k,])%*%invVMat%*%VB_MeanV[k,]-0.5*sum(diag(invVMat%*%VB_covMatVing))+psixi-sum(log(diag(VMat_chol)))
                                
                                ELBO2<- (obsMinHalf)%*%X_tilde%*%VB_MeanV[k,]+sum(log(diag(VB_covMatVchol))) #same as 1/2*log(det(VB_covMatVing))
                                
                                ELBO3<-VB_alpha_sigma*log(VB_beta_sigma)-log(gamma(VB_alpha_sigma))-VB_iter[k,1]*(log(VB_iter[k,2])) 
                                
                                ELBO4<- (VB_alpha_sigma-VB_iter[k,1])*(digamma(VB_iter[k,1])-log(VB_iter[k,2]))+(VB_iter[k,2]-VB_beta_sigma)*VB_E_sigma2Inv
                                
                                ELBOjthcal[k] <- ELBO1+ELBO2+ELBO3+ELBO4
                                
                                #Making ELBO for the Stopping Criteria
                                bar<-abs(ELBOjthcal[k]-ELBOjthcal[k-1])
                              }
                              return(c(ELBOjthcal[k],VB_MeanV[k,],VB_iter[k,],VB_covMatV[k,]))
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
thetaC_MeanV<-matrix(NA, nrow=iter.thetaphi, ncol=((0.8*n)+p))
dim(thetaC_MeanV)
thetaC_sigma2<-matrix(NA, nrow=iter.thetaphi, ncol=2)
dim(thetaC_sigma2)
thetaC_covMatV<-matrix(NA, nrow=iter.thetaphi, ncol=((0.8*n)+p))
dim(thetaC_covMatV)


for ( i in 1:iter.thetaphi){
  thetaC_MeanV[i,]<-outputMat[2:((0.8*n)+p+1),i]
}

for ( i in 1:iter.thetaphi){
  thetaC_sigma2[i,]<-outputMat[((0.8*n)+p+2):((0.8*n)+p+3),i]
}

for ( i in 1:iter.thetaphi){
  thetaC_covMatV[i,]<-outputMat[((0.8*n)+p+4):(1+((0.8*n)+p)+2+((0.8*n)+p)),i]
}


########################################################################
########################################################################
# Making new Density according to ELBOvector rates
########################################################################
########################################################################
# Beta1 Weighted Density
par(mfrow=c(1,1))
matDensityBeta1<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqb1<-seq(-3,3, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesBeta1 <-dnorm(x=xSeqb1, mean=thetaC_MeanV[k,1], sd=sqrt(thetaC_covMatV[k,1]))
  matDensityBeta1[,k]<-densitiesBeta1
}
dim(matDensityBeta1)
newdensityBeta1<-matDensityBeta1%*%ELBOvectorweights
#plot(x=xSeqb1, y=newdensityBeta1, typ="l", pch=16, col="red",main="Beta1")

# Beta2 Weighted Density
matDensityBeta2<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqb2<-seq(-3,3, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesBeta2 <-dnorm(x=xSeqb2, mean=thetaC_MeanV[k,p], sd=sqrt(thetaC_covMatV[k,2]))
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


############################################################
## W(Spatial random Effects)
############################################################
# W3 Weighted Density
q=3
matDensityW3<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqW3<-seq(-3,3, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesW3 <-dnorm(x=xSeqW3, mean=thetaC_MeanV[k,(p+q)], sd=sqrt(thetaC_covMatV[k,p+q]))
  matDensityW3[,k]<-densitiesW3
}
dim(matDensityW3)
newdensityW3<-matDensityW3%*%ELBOvectorweights


# W10 Weighted Density
q=10
matDensityW10<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqW10<-seq(-3,3, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesW10 <-dnorm(x=xSeqW10, mean=thetaC_MeanV[k,(p+q)], sd=sqrt(thetaC_covMatV[k,p+q]))
  matDensityW10[,k]<-densitiesW10
}
dim(matDensityW10)
newdensityW10<-matDensityW10%*%ELBOvectorweights


# W25 Weighted Density
q=25
matDensityW25<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqW25<-seq(-3,3, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesW25 <-dnorm(x=xSeqW25, mean=thetaC_MeanV[k,(p+q)], sd=sqrt(thetaC_covMatV[k,p+q]))
  matDensityW25[,k]<-densitiesW25
}
dim(matDensityW25)
newdensityW25<-matDensityW25%*%ELBOvectorweights


# W40 Weighted Density
q=40
matDensityW40<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
xSeqW40<-seq(-3,3, length.out=iter.thetaphi)
for ( k in 1:iter.thetaphi){
  densitiesW40 <-dnorm(x=xSeqW40, mean=thetaC_MeanV[k,(p+q)], sd=sqrt(thetaC_covMatV[k,p+q]))
  matDensityW40[,k]<-densitiesW40
}
dim(matDensityW40)
newdensityW40<-matDensityW40%*%ELBOvectorweights


# ########################################################################
# ########################################################################
# ## Prediction_VB
# ########################################################################
# ########################################################################
# VB_iter[j,1]
exp_prediction_VB<-matrix(rep(NA,(0.2*n)*iter.thetaphi), nrow=(0.2*n), ncol=iter.thetaphi)
eta_trainingdata_VB<-matrix(rep(NA,(0.8*n)*iter.thetaphi), nrow=(0.8*n), ncol=iter.thetaphi)
eta_testdata_VB<-matrix(rep(NA,(0.2*n)*iter.thetaphi), nrow=(0.2*n), ncol=iter.thetaphi)

for(i in 1:iter.thetaphi){
  if(i%%100==0) {print(i)}
  Sigma_11<-(thetaC_sigma2[i,2]/(thetaC_sigma2[i,1]-1))*exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/thetaphi[i])
  inv_Sigma_11<-solve(Sigma_11)
  Sigma_21<-(thetaC_sigma2[i,2]/(thetaC_sigma2[i,1]-1))*exp(-distMat[c((0.8*n+1):(n)),c(1:(0.8*n))]/thetaphi[i])
  Sigma_22<-(thetaC_sigma2[i,2]/(thetaC_sigma2[i,1]-1))*exp(-distMat[c((0.8*n+1):(n)),c((0.8*n+1):(n))]/thetaphi[i])
  Sigma_12<-(thetaC_sigma2[i,2]/(thetaC_sigma2[i,1]-1))*exp(-distMat[c(1:(0.8*n)),c((0.8*n+1):(n))]/thetaphi[i])
  eta_trainingdata_VB[,i]<-X_trainingdata%*%thetaC_MeanV[i,c(1,2)]+thetaC_MeanV[i,c((p+1):(0.8*n+p))]
  eta_testdata_VB[,i]<-X_testdata%*%thetaC_MeanV[i,c(1,2)]+Sigma_21%*%(inv_Sigma_11%*%thetaC_MeanV[i,c((p+1):(0.8*n+p))])
  exp_prediction_VB[,i]<-exp(eta_testdata_VB[,i])/(1+exp(eta_testdata_VB[,i]))
}

exp_prediction_VB<-exp_prediction_VB%*%ELBOvectorweights
obs_RMSPE_VB<-sqrt(mean((obs_testdata-exp_prediction_VB)^2))
obs_RMSPE_VB

#########################################################################
#########################################################################
obs_RMSPE_VB
VB_ptFinal
rocBinary <- roc(obs_testdata, exp_prediction_VB)
aucVal_VB<-auc(rocBinary)

save(exp_prediction_VB , obs_RMSPE_VB , aucVal_VB ,VB_ptFinal,newdensityBeta1, newdensityBeta2,
     newdensitySigma2,thetaphi,iter.thetaphi,ELBOvectorweights,
     newdensityW3, newdensityW10, newdensityW25, newdensityW40,
     thetaC_MeanV,thetaC_sigma2,thetaC_covMatV,
     file=paste0("B_phi",phiselect,"_Gen500_",sim,"_VBphiParallel.RData")) 
  }
}
