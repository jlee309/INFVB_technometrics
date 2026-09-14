################################################################################
################################################################################
# 10142024
# Running Multiple Simulations for Binary Data
# VB_Parallel : Fixing Phi and sigma2
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


## Preliminaries
# New functions for performing quadratic forms
################################################################################
quad.diag.New<-function (M, x) {
  colSums(crossprod(M, x) * x)
}
quad.form.new<-function(M,x){
  crossprod(crossprod(M, x), x)
}
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

## Initial Settings for VB
p<-length(beta)
iter<-1000 # (Note) Different from n.iteration in previous MCMC

## Beta
VB_mean_beta<-c(0,0)
VB_covMat_beta<-diag(2)
VB_inv_Sigma_beta_prior<-solve(VB_covMat_beta)

## Prior Distribution For Sigma2 (Inverse Gamma) 
VB_alpha_sigma<-0.1 ; VB_beta_sigma<-0.1
VB_E_sigma2Inv<-VB_alpha_sigma/VB_beta_sigma

## Discretizing ThetaD:{phi, sigma2} 
iter.thetaphissqphi=100 #Number of Discretizing Theta and sigma2
thetaphissqphi<-seq(0.001,sqrt(2),length.out=iter.thetaphissqphi)

x_sigma2=seq(0.00001,6,length.out=iter.thetaphissqphi)
theta_sigma2=dinvgamma(x=x_sigma2,shape=VB_alpha_sigma, rate=VB_beta_sigma,log=FALSE)
#plot(x_sigma2,theta_sigma2,col="blue")

VB_iteration<-iter.thetaphissqphi*iter.thetaphissqphi


################################################################################
################################################################################
## ELBOvector 
ELBOvector<-ELBOjthcal<-numeric()

## W 
R_phi_trainingdata<-exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/phi)
VB_invR_phi_trainingdata<-chol2inv(chol(R_phi_trainingdata))
dim(VB_invR_phi_trainingdata)

## Making X_tilde matrix 
X_tilde<-cbind(X_trainingdata,diag(0.8*n))
tX_tilde<-t(X_tilde)
dim(X_tilde)

## Making V(beta+W) Matrix
VB_MeanbetaMat<-matrix(NA, nrow=VB_iteration, ncol=p) 
dim(VB_MeanbetaMat)
VB_mean_wMat<-matrix(NA,nrow=VB_iteration,ncol=(0.8*n)) 
dim(VB_mean_wMat)

#First Line
VB_MeanbetaMat[1,]<-c(rep(1,p))
VB_mean_wMat[1,]<-c(rep(0,(0.8*n)))

#Making VB_MeanV Matrix
VB_MeanV<-cbind(VB_MeanbetaMat, VB_mean_wMat)
dim(VB_MeanV)

#Making covMatV Matrix 
VB_covMatV<-matrix(NA, nrow=VB_iteration, ncol=((0.8*n)+p))
dim(VB_covMatV)
VB_covMatV[1,]<-rep(1,((0.8*n)+p))
zeroMatpn<-matrix(0, nrow=p, ncol=(0.8*n)) 
zeroMatnp<-matrix(0, nrow=(0.8*n), ncol=p)
dim(zeroMatnp)

#Prior
invVMatmaking1<-cbind(VB_inv_Sigma_beta_prior,zeroMatpn)
dim(invVMatmaking1)
invVMatmaking2<-cbind(zeroMatnp,(1/sigma2)*VB_invR_phi_trainingdata)
dim(invVMatmaking2)
invVMat<-rbind(invVMatmaking1,invVMatmaking2)
dim(invVMat)


# Observations
obsMinHalf<-(obs_trainingdata-0.5)


## Final Saving Matrix for V(Beta, W) and sigma2
thetaC_MeanVssqphi<-matrix(NA, nrow=iter.thetaphissqphi^2, ncol=((0.8*n)+p))
dim(thetaC_MeanVssqphi)
thetaC_covMatVssqphi<-matrix(NA, nrow=iter.thetaphissqphi^2, ncol=(0.8*n+p))
dim(thetaC_covMatVssqphi)

#making temporary lambda to get xi 
foo<-(VB_covMatV[1,]+VB_MeanV[1,]%*%t(VB_MeanV[1,]))
xi<-sqrt(quad.diag.New(M=foo,x = tX_tilde))
lambda<- -tanh(xi/2)/(4*xi)
lambda<-diag(x = lambda)
tXtilLambdaXtil<-quad.form.new(M=lambda,x=X_tilde)


################################################################################
################################################################################
# INFVB using Parallelization Computing 
# Outer foreach / Inner for loop 
VB_pt<-proc.time()
outputMat<-foreach::foreach(j=1:iter.thetaphissqphi,.combine ="cbind",
                            .packages = c("mvtnorm","invgamma","fields")) %dopar% { 
                              
                              # Updating phi
                              R_phi_trainingdata<-exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/thetaphissqphi[j]) 
                              cholR_trainingdata<-chol(R_phi_trainingdata)
                              VB_invR_phi_trainingdata<-chol2inv(cholR_trainingdata) 
                              
                              for (k in 1:iter.thetaphissqphi){ 
                                # Updatng invVMat 
                                invVMatmaking2<-cbind(zeroMatnp,(1/x_sigma2[k])*VB_invR_phi_trainingdata)
                                invVMat<-rbind(invVMatmaking1,invVMatmaking2)
                                
                                VB_covMatchol<-chol(invVMat -2*tXtilLambdaXtil)
                                VB_covMatVing<-chol2inv(VB_covMatchol)
                                VB_covMatV[k,]<- diag(VB_covMatVing)
                                VB_MeanV[k,] <- VB_covMatVing%*%(tX_tilde%*%(obsMinHalf)) 
                                
                                # Update xi
                                foo<-(VB_covMatVing+VB_MeanV[k,]%*%t(VB_MeanV[k,]))
                                xi<-sqrt(quad.diag.New(M=foo,x = tX_tilde))
                                lambda<- -tanh(xi/2)/(4*xi)
                                lambda<-diag(x = lambda)
                                tXtilLambdaXtil<-quad.form.new(M=lambda,x=X_tilde)
                                
                                #Updata psi(xi)
                                psixi<-xi/2-log(1+exp(xi))+xi*tanh(xi/2)/4
                                psixi<-sum(psixi)
                                
                                ################################################################################
                                # ELBO preparation
                                VB_covMatVchol<-chol(VB_covMatVing)
                                invVMat_chol<-chol(invVMat)
                                VMat<-chol2inv(invVMat_chol)
                                VMat_chol<-chol(VMat)
                                
                                # ELBO Calculation
                                ELBO1<-t(VB_MeanV[k,])%*%tX_tilde%*%lambda%*%(X_tilde)%*%VB_MeanV[k,]+sum(diag(tX_tilde%*%(lambda%*%(X_tilde))%*%VB_covMatVing))-0.5*t(VB_MeanV[k,])%*%invVMat%*%VB_MeanV[k,]-0.5*sum(diag(invVMat%*%VB_covMatVing))
                                
                                ELBO2<- (obsMinHalf)%*%X_tilde%*%VB_MeanV[k,]+psixi-sum(log(diag(VMat_chol)))
                                
                                ELBO3<- sum(log(diag(VB_covMatVchol))) +(0.8*n+p)/2 #same as 1/2*log(det(VB_covMatVing))
                                
                                ELBO4<- log(theta_sigma2[k])
                                
                                ELBOjthcal[k] <- ELBO1+ELBO2+ELBO3+ELBO4
                              }
                              return(c(ELBOjthcal[1:k],VB_MeanV[(1:k),],VB_covMatV[(1:k),]))
                            }
VB_ptFinalssqphi<-proc.time()-VB_pt

################################################################################
################################################################################
#Getting ELBO 
ELBOvector<-ELBOvectorcal<-numeric()
for(k in 1:iter.thetaphissqphi){
  ELBOvectorcal<-outputMat[1:iter.thetaphissqphi,k]
  ELBOvector<-c(ELBOvector,ELBOvectorcal)
}
weights<-exp(ELBOvector-max(ELBOvector))
ELBOvectorweights<-weights/sum(weights)


# Making tables for ELBOs  
thetamatrix<-expand.grid(x_sigma2,thetaphissqphi) 
head(thetamatrix,25)
dim(thetamatrix)
ELBOtables<-cbind(thetamatrix,ELBOvectorweights)
dim(ELBOtables)

# ELBO for phi
ELBOphi<-c()
for(i in 1:iter.thetaphissqphi){
  ELBOphiprev<-sum(subset(ELBOtables, Var2==thetaphissqphi[i])$ELBOvectorweights)
  ELBOphi<-c(ELBOphi,ELBOphiprev)
}

# ELBO for sigma2
ELBOsigma2<-c()
for(i in 1:iter.thetaphissqphi){
  ELBOsigma2prev<-sum(subset(ELBOtables, Var1==x_sigma2[i])$ELBOvectorweights)
  ELBOsigma2<-c(ELBOsigma2,ELBOsigma2prev)
}

################################################################################
################################################################################
#Getting thetaC_MeanVssqphi
dim(thetaC_MeanVssqphi)
for(k in 1:iter.thetaphissqphi){
  for(i in 1:iter.thetaphissqphi){
    rows<-numeric()
    for(number in 1:dim(VB_MeanV)[2]){
      rows[number]<-(i+iter.thetaphissqphi*(number))
    }
    thetaC_MeanVssqphi[i+iter.thetaphissqphi*(k-1),]<-outputMat[rows,k]
  }
}


################################################################################
################################################################################
#Getting thetaC_covMatVssqphi
dim(thetaC_covMatVssqphi)

for(k in 1:iter.thetaphissqphi){
  for(i in 1:iter.thetaphissqphi){
    rows<-numeric()
    for(number in 1:dim(VB_MeanV)[2]){
      rows[number]<-(i+(0.8*n+p)*iter.thetaphissqphi+iter.thetaphissqphi*(number))
    }
    thetaC_covMatVssqphi[i+iter.thetaphissqphi*(k-1),]<-outputMat[rows,k]
  }
}

################################################################################
################################################################################
################################################################################
################################################################################
# Making new Density according to ELBOvector rates
################################################################################
################################################################################
# Beta1 Weighted Density 
par(mfrow=c(1,1))
matDensityBeta1<-matrix(NA,ncol=iter.thetaphissqphi^2,nrow=iter.thetaphissqphi^2)
xSeqb1<-seq(-3,3, length.out=iter.thetaphissqphi^2)
for ( k in 1:iter.thetaphissqphi^2){
  densitiesBeta1 <-dnorm(x=xSeqb1, mean=thetaC_MeanVssqphi[k,1], sd=sqrt(thetaC_covMatVssqphi[k,1]))
  matDensityBeta1[,k]<-densitiesBeta1
}
dim(matDensityBeta1)
newdensityBeta1ssqphi<-matDensityBeta1%*%ELBOvectorweights
#plot(x=xSeqb1, y=newdensityBeta1ssqphi, typ="l", pch=16, col="red",main="Beta1")



# Beta2 Weighted Density 
matDensityBeta2<-matrix(NA,ncol=iter.thetaphissqphi^2,nrow=iter.thetaphissqphi^2)
xSeqb2<-seq(-3,3, length.out=iter.thetaphissqphi^2)
for ( k in 1:iter.thetaphissqphi^2){
  densitiesBeta2 <-dnorm(x=xSeqb2, mean=thetaC_MeanVssqphi[k,p], sd=sqrt(thetaC_covMatVssqphi[k,2]))
  matDensityBeta2[,k]<-densitiesBeta2
}
dim(matDensityBeta2)
newdensityBeta2ssqphi<-matDensityBeta2%*%ELBOvectorweights
#plot(x=xSeqb2, y=newdensityBeta2ssqphi, typ="l", pch=16, col="red", xlim=c(-3,3), main="Beta2")

################################################################################
################################################################################
################################################################################
################################################################################
## W(Spatial random Effects)
################################################################################
# W3 Weighted Density (Wqth) 
q=3
matDensityW3<-matrix(NA,ncol=iter.thetaphissqphi^2,nrow=iter.thetaphissqphi^2)
xSeqW3<-seq(-3,3, length.out=iter.thetaphissqphi^2)
for ( k in 1:iter.thetaphissqphi^2){
  densitiesW3 <-dnorm(x=xSeqW3, mean=thetaC_MeanVssqphi[k,(p+q)], sd=sqrt(thetaC_covMatVssqphi[k,(p+q)]))
  matDensityW3[,k]<-densitiesW3
}
dim(matDensityW3)
newdensityW3ssqphi<-matDensityW3%*%ELBOvectorweights
#plot(x=xSeqW3, y=newdensityW3ssqphi, typ="l", pch=16, col="red", xlim=c(-3,3), main="")

# W10 Weighted Density (Wqth) 
q=10
matDensityW10<-matrix(NA,ncol=iter.thetaphissqphi^2,nrow=iter.thetaphissqphi^2)
xSeqW10<-seq(-3,3, length.out=iter.thetaphissqphi^2)
for ( k in 1:iter.thetaphissqphi^2){
  densitiesW10 <-dnorm(x=xSeqW10, mean=thetaC_MeanVssqphi[k,(p+q)], sd=sqrt(thetaC_covMatVssqphi[k,(p+q)]))
  matDensityW10[,k]<-densitiesW10
}
dim(matDensityW10)
newdensityW10ssqphi<-matDensityW10%*%ELBOvectorweights

# W25 Weighted Density (Wqth) 
q=25
matDensityW25<-matrix(NA,ncol=iter.thetaphissqphi^2,nrow=iter.thetaphissqphi^2)
xSeqW25<-seq(-3,3, length.out=iter.thetaphissqphi^2)
for ( k in 1:iter.thetaphissqphi^2){
  densitiesW25 <-dnorm(x=xSeqW25, mean=thetaC_MeanVssqphi[k,(p+q)], sd=sqrt(thetaC_covMatVssqphi[k,(p+q)]))
  matDensityW25[,k]<-densitiesW25
}
dim(matDensityW25)
newdensityW25ssqphi<-matDensityW25%*%ELBOvectorweights

# W40 Weighted Density (Wqth) 
q=40
matDensityW40<-matrix(NA,ncol=iter.thetaphissqphi^2,nrow=iter.thetaphissqphi^2)
xSeqW40<-seq(-3,3, length.out=iter.thetaphissqphi^2)
for ( k in 1:iter.thetaphissqphi^2){
  densitiesW40 <-dnorm(x=xSeqW40, mean=thetaC_MeanVssqphi[k,(p+q)], sd=sqrt(thetaC_covMatVssqphi[k,(p+q)]))
  matDensityW40[,k]<-densitiesW40
}
dim(matDensityW40)
newdensityW40ssqphi<-matDensityW40%*%ELBOvectorweights


# ################################################################################
# ################################################################################
# # Eta1 Finding ---> We will do a sampling method
# ################################################################################
# ################################################################################
# # Eta Mean Vector
# EXeta_Meanssqphi<-(X_tilde%*%t(thetaC_MeanVssqphi))%*%ELBOvectorweights
# dim(EXeta_Meanssqphi)
# 
# ################################################################################
# # Eta Variance-Covariance Matrix 
# q=3
# etai_CovMatssqphi<-matrix(data=NA, nrow=1, ncol=iter.thetaphissqphi^2)
# for(i in 1:iter.thetaphissqphi^2){
#   etamat<-matrix(data=thetaC_covMatVssqphi[i,], nrow=(n*0.8)+p, ncol=(n*0.8)+p, byrow = TRUE)  
#   Upperetamat<-cbind(etamat[(1:2),(1:2)],etamat[(1:2),(p+q)])
#   Loweretamat<-c(etamat[(p+q),(1:2)],etamat[(p+q),(p+q)])
#   etamatcal<-rbind(Upperetamat,Loweretamat)
#   etamatcal<-as.matrix(etamatcal)
#   X_tilde_q<-cbind(t(X_trainingdata[q,]), diag(1))
#   etai_CovMatssqphi[,i]<-X_tilde_q%*%etamatcal%*%t(X_tilde_q)
# }
# EXetai_CovMatssqphi_3<-etai_CovMatssqphi%*%ELBOvectorweights
# 
# ## Get only the highest value
# ELBOvectorweightsSq<-ifelse(ELBOvectorweights==max(ELBOvectorweights),1,0 )
# # ELBOvectorweightsSq<-ELBOvectorweights^2
# EXetai_CovMatssqphi_3_Sq<-etai_CovMatssqphi%*%ELBOvectorweightsSq
# 
# 
# ## Get the two highest values
# ELBOvectorweightsSq_toptwo <- sort(ELBOvectorweights, decreasing = TRUE)[1:2]
# # Replace all values except the two highest with 0
# rescaled_vector <- ifelse(ELBOvectorweights %in% ELBOvectorweightsSq_toptwo, ELBOvectorweights, 0)
# # Rescale the vector so the sum of the elements equals 1
# ELBOvectorweightsSq_toptwo <- rescaled_vector / sum(rescaled_vector)
# 
# EXetai_CovMat_3_Sq_toptwo_ssqphi<-etai_CovMatssqphi%*%ELBOvectorweightsSq_toptwo
# 
# 
# 
# 
# # Eta Variance-Covariance Matrix 
# q=10
# etai_CovMatssqphi<-matrix(data=NA, nrow=1, ncol=iter.thetaphissqphi^2)
# for(i in 1:iter.thetaphissqphi^2){
#   etamat<-matrix(data=thetaC_covMatVssqphi[i,], nrow=(n*0.8)+p, ncol=(n*0.8)+p, byrow = TRUE)  
#   Upperetamat<-cbind(etamat[(1:2),(1:2)],etamat[(1:2),(p+q)])
#   Loweretamat<-c(etamat[(p+q),(1:2)],etamat[(p+q),(p+q)])
#   etamatcal<-rbind(Upperetamat,Loweretamat)
#   etamatcal<-as.matrix(etamatcal)
#   X_tilde_q<-cbind(t(X_trainingdata[q,]), diag(1))
#   etai_CovMatssqphi[,i]<-X_tilde_q%*%etamatcal%*%t(X_tilde_q)
# }
# EXetai_CovMatssqphi_10<-etai_CovMatssqphi%*%ELBOvectorweights
# EXetai_CovMatssqphi_10_Sq<-etai_CovMatssqphi%*%ELBOvectorweightsSq
# EXetai_CovMat_10_Sq_toptwo_ssqphi<-etai_CovMatssqphi%*%ELBOvectorweightsSq_toptwo
# 
# # Eta Variance-Covariance Matrix 
# q=25
# etai_CovMatssqphi<-matrix(data=NA, nrow=1, ncol=iter.thetaphissqphi^2)
# for(i in 1:iter.thetaphissqphi^2){
#   etamat<-matrix(data=thetaC_covMatVssqphi[i,], nrow=(n*0.8)+p, ncol=(n*0.8)+p, byrow = TRUE)  
#   Upperetamat<-cbind(etamat[(1:2),(1:2)],etamat[(1:2),(p+q)])
#   Loweretamat<-c(etamat[(p+q),(1:2)],etamat[(p+q),(p+q)])
#   etamatcal<-rbind(Upperetamat,Loweretamat)
#   etamatcal<-as.matrix(etamatcal)
#   X_tilde_q<-cbind(t(X_trainingdata[q,]), diag(1))
#   etai_CovMatssqphi[,i]<-X_tilde_q%*%etamatcal%*%t(X_tilde_q)
# }
# EXetai_CovMatssqphi_25<-etai_CovMatssqphi%*%ELBOvectorweights
# EXetai_CovMatssqphi_25_Sq<-etai_CovMatssqphi%*%ELBOvectorweightsSq
# EXetai_CovMat_25_Sq_toptwo_ssqphi<-etai_CovMatssqphi%*%ELBOvectorweightsSq_toptwo
# 
# 
# # Eta Variance-Covariance Matrix 
# q=40
# etai_CovMatssqphi<-matrix(data=NA, nrow=1, ncol=iter.thetaphissqphi^2)
# for(i in 1:iter.thetaphissqphi^2){
#   etamat<-matrix(data=thetaC_covMatVssqphi[i,], nrow=(n*0.8)+p, ncol=(n*0.8)+p, byrow = TRUE)  
#   Upperetamat<-cbind(etamat[(1:2),(1:2)],etamat[(1:2),(p+q)])
#   Loweretamat<-c(etamat[(p+q),(1:2)],etamat[(p+q),(p+q)])
#   etamatcal<-rbind(Upperetamat,Loweretamat)
#   etamatcal<-as.matrix(etamatcal)
#   X_tilde_q<-cbind(t(X_trainingdata[q,]), diag(1))
#   etai_CovMatssqphi[,i]<-X_tilde_q%*%etamatcal%*%t(X_tilde_q)
# }
# EXetai_CovMatssqphi_40<-etai_CovMatssqphi%*%ELBOvectorweights
# EXetai_CovMatssqphi_40_Sq<-etai_CovMatssqphi%*%ELBOvectorweightsSq
# EXetai_CovMat_40_Sq_toptwo_ssqphi<-etai_CovMatssqphi%*%ELBOvectorweightsSq_toptwo
# 


################################################################################
################################################################################
## Prediction_VB
################################################################################
################################################################################
exp_prediction_VBssqphi<-matrix(rep(NA,(0.2*n)*iter.thetaphissqphi^2), nrow=(0.2*n), ncol=iter.thetaphissqphi^2)
eta_trainingdata_VBssqphi<-matrix(rep(NA,(0.8*n)*iter.thetaphissqphi^2), nrow=(0.8*n), ncol=iter.thetaphissqphi^2)
eta_testdata_VBssqphi<-matrix(rep(NA,(0.2*n)*iter.thetaphissqphi^2), nrow=(0.2*n), ncol=iter.thetaphissqphi^2)


for(i in 1:iter.thetaphissqphi^2){
  if(i%%100==0) {print(i)}
  Sigma_11<-(thetamatrix[i,1])*exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/thetamatrix[i,2])
  inv_Sigma_11<-solve(Sigma_11)
  Sigma_21<-(thetamatrix[i,1])*exp(-distMat[c((0.8*n+1):(n)),c(1:(0.8*n))]/thetamatrix[i,2])
  eta_trainingdata_VBssqphi[,i]<-X_trainingdata%*%thetaC_MeanVssqphi[i,c(1,2)]+thetaC_MeanVssqphi[i,c((p+1):(0.8*n+p))]
  eta_testdata_VBssqphi[,i]<-X_testdata%*%thetaC_MeanVssqphi[i,c(1,2)]+Sigma_21%*%(inv_Sigma_11%*%thetaC_MeanVssqphi[i,c((p+1):(0.8*n+p))])
  exp_prediction_VBssqphi[,i]<-exp(eta_testdata_VBssqphi[,i])/(1+exp(eta_testdata_VBssqphi[,i]))
}

dim(exp_prediction_VBssqphi)
exp_prediction_VBssqphi<-exp_prediction_VBssqphi%*%ELBOvectorweights
obs_RMSPE_VBssqphi<-sqrt(mean((obs_testdata-exp_prediction_VBssqphi)^2))
obs_RMSPE_VBssqphi

################################################################################
################################################################################
rocBinary <- roc(obs_testdata, exp_prediction_VBssqphi)
aucVal_VBssqphi<-auc(rocBinary)

save(exp_prediction_VBssqphi , obs_RMSPE_VBssqphi , aucVal_VBssqphi ,VB_ptFinalssqphi,
     ELBOsigma2,ELBOphi,
     newdensityBeta1ssqphi, newdensityBeta2ssqphi,x_sigma2,thetaphissqphi,iter.thetaphissqphi,
     newdensityW3ssqphi,newdensityW10ssqphi,newdensityW25ssqphi,newdensityW40ssqphi,ELBOvectorweights,
     thetaC_MeanVssqphi,thetaC_covMatVssqphi,ELBOtables,
    file = paste0("B_phi",phiselect,"_Gen500_",sim,"_VBphissqphiParallel.RData"))
      }
    }

