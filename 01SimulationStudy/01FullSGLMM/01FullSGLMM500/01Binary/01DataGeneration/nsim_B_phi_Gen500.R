################################################################################
################################################################################
# 10142024
# Binary Data Generation when N=500
# Estimating V(Beta,W),Sigma2,phi 
################################################################################
################################################################################
getwd()
rm(list=ls())
#dev.off()
library(cli);library(viridis);library(fields);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm);
library(invgamma);library(pROC);library(Matrix);library(emulator);
################################################################################
################################################################################
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
# set.seed(123)  # Set an initial seed
set.seed(5153)
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed

for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    n=500 # Sample Size
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    ## Parameters
    beta<-c(1,1) ; sigma2=1 ; phi=phiset[phiselect] ; p<-length(beta)
    ## Locations and Design Matrix (X)
    X<-cbind(runif(n,-1,1), runif(n,-1,1))
    locations<-cbind(runif(n),runif(n))
    distMat<-as.matrix(rdist(locations)) 
    R_phi<-sigma2*exp(-distMat/phi) 
    w<-t(chol(R_phi))%*%rnorm(n) # Spatial Random Effects W
    XB<-X%*%beta
    eta<-XB+w
    prob<-exp(eta)/(1+exp(eta))
    summary(prob)
    obs<-rbinom(n=n, size=1, prob=prob)
    
    ################################################################################
    ################################################################################
    ## Splitting into Two Data(Training, Test Data)
    ################################################################################
    ################################################################################
    #Training Data (80% of the whole data)
    X_trainingdata<-X[c(1:(0.8*n)),]
    locations_trainingdata<-locations[c(1:(0.8*n)),]
    obs_trainingdata<-obs[c(1:(0.8*n))]
    w_trainingdata<-w[c(1:((0.8*n))),]
    prob_trainingdata<-prob[c(1:((0.8*n))),]
    R_phi_trainingdata<-exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/phi)
    sigma2*R_phi_trainingdata[1,40]==R_phi[1,40] #Checking the R_phi_trainingdata and the full R_phi
     
    # #Test Data (20% of the whole data)
    X_testdata<-X[c(((0.8*n)+1):(n)),]
    locations_testdata<-locations[c(((0.8*n)+1):(n)),]
    obs_testdata<-obs[c(((0.8*n)+1):(n))]
    w_testdata<-w[c((0.8*n+1):(n)),]
    prob_testdata<-prob[c(((0.8*n)+1):(n)),]
    ################################################################################
    ################################################################################
    save(phi,prob,p,w,distMat,beta, sigma2, phi, n,
              X_trainingdata,locations_trainingdata,obs_trainingdata,w_trainingdata,
              prob_trainingdata,R_phi_trainingdata,
              X_testdata, locations_testdata, obs_testdata, w_testdata,prob_testdata,
         file=paste0("B_phi",phiselect,"_Gen500_",sim,".RData")) 
  } 
}

################################################################################
################################################################################


