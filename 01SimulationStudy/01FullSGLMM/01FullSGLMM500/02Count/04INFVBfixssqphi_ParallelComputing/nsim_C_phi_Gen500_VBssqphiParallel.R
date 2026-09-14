################################################################################
################################################################################
# 10142024
# Running Multiple Simulations for Count Data
# VB_Parallel : Fixing Phi and sigma2
################################################################################
################################################################################


################################################################################
################################################################################
# Variational Bayes (INFVB) Fixing Sigma2 and phi, Count data when phi=0.1
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
## INFVB 
# thetaC={Gamma(Beta,W)} : Estimating Gamma, Sigma2 and Phi
# thetaD={phi,Sigma2} : Discretizing this Phi and Sigma2; 
################################################################################
################################################################################
n_simulations <- 50  #Number of simulations per each scenario
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


## Preliminaries
p<-length(beta)
iter<-1000 # (Note) Different from n.iteration in previous MCMC
## VB Beta Prior (Same as MCMC prior)
VB_mu_beta_prior<-c(0,0); VB_sigma_beta_prior<-diag(2) 
chol_VB_sigma_beta_prior<-chol(VB_sigma_beta_prior)
VB_inv_Sigma_beta_prior<-chol2inv(chol_VB_sigma_beta_prior)

## VB Sigma2 Prior
#VB_iter<-matrix(NA, nrow=iter, ncol=2) 
#VB_iter[1,]<-c(1,1) #Initial values
VB_alpha_sigma<-2 ; VB_beta_sigma<-2
#VB_E_sigma2Inv<-VB_iter[1,1]/VB_iter[1,2]

## Discretizing ThetaD:{phi, sigma2} 
iter.thetaphissqphi=100 #Number of Discretizing Theta and sigma2
thetaphissqphi<-seq(0.001,sqrt(2),length.out=iter.thetaphissqphi)

x_sigma2=seq(0.00001,6,length.out=iter.thetaphissqphi)
theta_sigma2=dinvgamma(x=x_sigma2,shape=VB_alpha_sigma, rate=VB_beta_sigma,log=FALSE)
#plot(x_sigma2,theta_sigma2,col="blue")

VB_iteration<-iter.thetaphissqphi^2

## ELBOvector 
ELBOvector<-ELBOjthcal<-numeric()

## Final Saving Matrix for Gamma(Beta, W)
thetaC_Meangammassqphi<-matrix(NA, nrow=iter.thetaphissqphi^2, ncol=((0.8*n)+p))
dim(thetaC_Meangammassqphi)
thetaC_covMatgammassqphi<-matrix(NA, nrow=iter.thetaphissqphi^2, ncol=(0.8*n+p))
dim(thetaC_covMatgammassqphi)

## VB_invR_phi for training data
R_phi_trainingdata<-exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/phi)
VB_invR_phi_trainingdata<-chol2inv(chol(R_phi_trainingdata))

## Making X_tilde(X+Identity)
X_tilde_trainingdata<-cbind(X_trainingdata,diag(0.8*n))
dim(X_tilde_trainingdata)

## Making Gamma(beta+W) Matrix
VB_MeanbetaMat<-matrix(NA, nrow=VB_iteration, ncol=2) 
VB_MeanWMat<-matrix(NA, nrow=VB_iteration, ncol=(0.8*n)) #Only for training W
VB_MeangammaS<-cbind(VB_MeanbetaMat, VB_MeanWMat)
VB_MeangammaS[1,]<-rep(1,(0.8*n+p))
dim(VB_MeangammaS)
VB_covMatgammaS<-matrix(NA, nrow=VB_iteration, ncol=(0.8*n+p))
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
# INFVB using Parallelization Computing 
# Outer foreach / Inner for loop 
VB_pt<-proc.time()
outputMat<-foreach::foreach(j=1:iter.thetaphissqphi,.combine ="cbind",
                            .packages = c("mvtnorm","invgamma","fields")) %dopar% { 
                              # Updating phi
                              R_phi_trainingdata<-exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/thetaphissqphi[j]) 
                              cholR_trainingdata<-chol(R_phi_trainingdata)
                              VB_invR_phi_trainingdata<-chol2inv(cholR_trainingdata) 
                              
                              
                              ELBOjthcal<-numeric()
                              ELBOjthcal[1]<-10
                              
                              
                              for(k in 1:iter.thetaphissqphi) {
                                
                                # Updatng Discreticized Phi and E[1/(sigma^2)] for invgammaMat and gammaMat
                                invgammaMatmaking2<-cbind(zeroMatn2,1/x_sigma2[k]*VB_invR_phi_trainingdata)
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
                                
                                #ELBO preparation
                                VB_covMatgamma_chol<-chol(VB_covMatgamma)
                                
                                # ELBO 
                                ELBO1<- t(obs_trainingdata)%*%(X_tilde_trainingdata%*%VB_Meangamma)-t(rep(1,(0.8*n)))%*%exp(X_tilde_trainingdata%*%VB_Meangamma+(1/2)*diag(X_tilde_trainingdata%*%VB_covMatgamma%*%t(X_tilde_trainingdata)))
                                
                                ELBO2<- -sum(log(diag(cholgammaMat)))-(1/2)*t(VB_Meangamma)%*%(invgammaMat)%*%(VB_Meangamma)-(1/2)*sum(diag(invgammaMat%*%VB_covMatgamma))
                                
                                ELBO3<- sum(log(diag(VB_covMatgamma_chol)))
                                
                                ELBO4<- log(theta_sigma2[k])
                                
                                ELBOjthcal[k] <- ELBO1+ELBO2+ELBO3+ELBO4
                              }
                              return(c(ELBOjthcal[1:k],VB_MeangammaS[(1:k),],VB_covMatgammaS[(1:k),]))
                            }
VB_ptFinalssqphi<-proc.time()-VB_pt
VB_ptFinalssqphi

dim(outputMat)

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
#Getting thetaC_Meangammassqphi
dim(thetaC_Meangammassqphi)
for(k in 1:iter.thetaphissqphi){
  for(i in 1:iter.thetaphissqphi){
    rows<-numeric()
    for(number in 1:dim(VB_MeangammaS)[2]){
      rows[number]<-(i+iter.thetaphissqphi*(number))
    }
    thetaC_Meangammassqphi[i+iter.thetaphissqphi*(k-1),]<-outputMat[rows,k]
  }
}


################################################################################
################################################################################
#Getting thetaC_covMatgammassqphi
dim(thetaC_covMatgammassqphi)
for(k in 1:iter.thetaphissqphi){
  for(i in 1:iter.thetaphissqphi){
    rows<-numeric()
    for(number in 1:dim(VB_MeangammaS)[2]){
      rows[number]<-(i+(0.8*n+p)*iter.thetaphissqphi+iter.thetaphissqphi*(number))
    }
    thetaC_covMatgammassqphi[i+iter.thetaphissqphi*(k-1),]<-outputMat[rows,k]
  }
}


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
  densitiesBeta1 <-dnorm(x=xSeqb1, mean=thetaC_Meangammassqphi[k,1], sd=sqrt(thetaC_covMatgammassqphi[k,1]))
  matDensityBeta1[,k]<-densitiesBeta1
}
dim(matDensityBeta1)
newdensityBeta1ssqphi<-matDensityBeta1%*%ELBOvectorweights
#plot(x=xSeqb1, y=newdensityBeta1ssqphi, typ="l", pch=16, col="red",main="Beta1")

# Beta2 Weighted Density 
matDensityBeta2<-matrix(NA,ncol=iter.thetaphissqphi^2,nrow=iter.thetaphissqphi^2)
xSeqb2<-seq(-3,3, length.out=iter.thetaphissqphi^2)
for ( k in 1:iter.thetaphissqphi^2){
  densitiesBeta2 <-dnorm(x=xSeqb2, mean=thetaC_Meangammassqphi[k,p], sd=sqrt(thetaC_covMatgammassqphi[k,2]))
  matDensityBeta2[,k]<-densitiesBeta2
}
dim(matDensityBeta2)
newdensityBeta2ssqphi<-matDensityBeta2%*%ELBOvectorweights
dim(newdensityBeta2ssqphi)
#plot(x=xSeqb2, y=newdensityBeta1ssqphi, typ="l", pch=16, col="red", xlim=c(0,2), main="Beta2")



################################################################################
################################################################################
## W(Spatial random Effects)
################################################################################
################################################################################
# W3 Weighted Density (Wqth) 
q=3
matDensityW3<-matrix(NA,ncol=iter.thetaphissqphi^2,nrow=iter.thetaphissqphi^2)
xSeqW3<-seq(-3,3, length.out=iter.thetaphissqphi^2)
for ( k in 1:iter.thetaphissqphi^2){
  densitiesW3 <-dnorm(x=xSeqW3, mean=thetaC_Meangammassqphi[k,(p+q)], sd=sqrt(thetaC_covMatgammassqphi[k,p+q]))
  matDensityW3[,k]<-densitiesW3
}
dim(matDensityW3)
newdensityW3ssqphi<-matDensityW3%*%ELBOvectorweights

# W10 Weighted Density (Wqth) 
q=10
matDensityW10<-matrix(NA,ncol=iter.thetaphissqphi^2,nrow=iter.thetaphissqphi^2)
xSeqW10<-seq(-3,3, length.out=iter.thetaphissqphi^2)
for ( k in 1:iter.thetaphissqphi^2){
  densitiesW10 <-dnorm(x=xSeqW10, mean=thetaC_Meangammassqphi[k,(p+q)], sd=sqrt(thetaC_covMatgammassqphi[k,p+q]))
  matDensityW10[,k]<-densitiesW10
}
dim(matDensityW10)
newdensityW10ssqphi<-matDensityW10%*%ELBOvectorweights

# W25 Weighted Density (Wqth) 
q=25
matDensityW25<-matrix(NA,ncol=iter.thetaphissqphi^2,nrow=iter.thetaphissqphi^2)
xSeqW25<-seq(-3,3, length.out=iter.thetaphissqphi^2)
for ( k in 1:iter.thetaphissqphi^2){
  densitiesW25 <-dnorm(x=xSeqW25, mean=thetaC_Meangammassqphi[k,(p+q)], sd=sqrt(thetaC_covMatgammassqphi[k,p+q]))
  matDensityW25[,k]<-densitiesW25
}
dim(matDensityW25)
newdensityW25ssqphi<-matDensityW25%*%ELBOvectorweights

# W40 Weighted Density (Wqth) 
q=40
matDensityW40<-matrix(NA,ncol=iter.thetaphissqphi^2,nrow=iter.thetaphissqphi^2)
xSeqW40<-seq(-3,3, length.out=iter.thetaphissqphi^2)
for ( k in 1:iter.thetaphissqphi^2){
  densitiesW40 <-dnorm(x=xSeqW40, mean=thetaC_Meangammassqphi[k,(p+q)], sd=sqrt(thetaC_covMatgammassqphi[k,p+q]))
  matDensityW40[,k]<-densitiesW40
}
dim(matDensityW40)
newdensityW40ssqphi<-matDensityW40%*%ELBOvectorweights




################################################################################
################################################################################
## Prediction_VB
################################################################################
################################################################################

lambda_prediction_VB<-matrix(rep(NA,(0.2*n)*iter.thetaphissqphi^2), nrow=(0.2*n), ncol=iter.thetaphissqphi^2)
eta_trainingdata_VBssqphi<-matrix(rep(NA,(0.8*n)*iter.thetaphissqphi^2), nrow=(0.8*n), ncol=iter.thetaphissqphi^2)
eta_testdata_VBssqphi<-matrix(rep(NA,(0.2*n)*iter.thetaphissqphi^2), nrow=(0.2*n), ncol=iter.thetaphissqphi^2)


for(i in 1:iter.thetaphissqphi^2){
  if(i%%100==0) {print(i)}
  Sigma_11<-(thetamatrix[i,1])*exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/thetamatrix[i,2])
  inv_Sigma_11<-solve(Sigma_11)
  Sigma_21<-(thetamatrix[i,1])*exp(-distMat[c((0.8*n+1):(n)),c(1:(0.8*n))]/thetamatrix[i,2])
  
  eta_trainingdata_VBssqphi[,i]<-X_trainingdata%*%thetaC_Meangammassqphi[i,c(1,2)]+thetaC_Meangammassqphi[i,c((p+1):(0.8*n+p))]
  eta_testdata_VBssqphi[,i]<-X_testdata%*%thetaC_Meangammassqphi[i,c(1,2)]+Sigma_21%*%(inv_Sigma_11%*%thetaC_Meangammassqphi[i,c((p+1):(0.8*n+p))])
  lambda_prediction_VB[,i]<-exp(eta_testdata_VBssqphi[,i])
}

lambda_prediction_VBssqphi<-lambda_prediction_VB%*%ELBOvectorweights
lambda_RMSPE_VBssqphi<-sqrt(mean((lambda_testdata-lambda_prediction_VBssqphi)^2))
obs_RMSPE_VBssqphi<-sqrt(mean((obs_testdata-lambda_prediction_VBssqphi)^2))



save(obs_RMSPE_VBssqphi,lambda_RMSPE_VBssqphi,VB_ptFinalssqphi,ELBOsigma2,ELBOphi,
     newdensityBeta1ssqphi, newdensityBeta2ssqphi,x_sigma2,thetaphissqphi,iter.thetaphissqphi,
     newdensityW3ssqphi,newdensityW10ssqphi,newdensityW25ssqphi,newdensityW40ssqphi,ELBOvectorweights,
     thetaC_Meangammassqphi,thetaC_covMatgammassqphi,ELBOtables,
     file = paste0("C_phi",phiselect,"_Gen500_",sim,"_VBphissqphiParallel.RData"))
  }
}

