################################################################################
################################################################################
# 10192024 
# Gaussian Data Generation When N=500 
# Estimating Beta,Sigma2,phi, not W 
################################################################################
################################################################################
getwd()
rm(list=ls())
#dev.off()
library(cli);library(viridis);library(fields);library(pgdraw);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm)
library(invgamma);library(pROC);library(Matrix);library(emulator)

################################################################################
################################################################################
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(5153)
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed

for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    n=500 # Number of Locations
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)

## Parameters
beta<-c(1,1) ; sigma2=1 ; phi=phiset[phiselect] # when phi increases, MCMC and VB looks much similar
p<-length(beta)
## Locations and Design Matrix(X)
X<-cbind(runif(n,-1,1), runif(n,-1,1)) # X are fixed
locations<-cbind(runif(n),runif(n))    # Locations are fixed
distMat<-as.matrix(rdist(locations))  # Distance Matrix
dim(distMat)
R_phi<-exp(-distMat/phi)
dim(R_phi)
w<-t(chol(R_phi))%*%rnorm(n, mean=0, sd=sqrt(sigma2)) #Spatial Random Effects W
XB<-X%*%beta
Z<-XB+w

################################################################################
################################################################################
## Splitting into Two Data(Training, Test Data)
################################################################################
################################################################################
#Training Data (80% of the whole data)
X_trainingdata<-X[c(1:(0.8*n)),]
locations_trainingdata<-locations[c(1:(0.8*n)),]
Z_trainingdata<-Z[c(1:(0.8*n)),]
w_trainingdata<-w[c(1:(0.8*n)),]
R_phi_trainingdata<-exp(-distMat[c(1:(0.8*n)),c(1:(0.8*n))]/phi)
dim(R_phi_trainingdata)
R_phi_trainingdata[1,40];R_phi[1,40] # Just want to check 


#Test Data (20% of the whole data)
X_testdata<-X[c(((0.8*n)+1):(n)),]
locations_testdata<-locations[c(((0.8*n)+1):(n)),]
Z_testdata<-Z[c(((0.8*n)+1):(n)),]
w_testdata<-w[c(((0.8*n)+1):(n)),]

################################################################################
################################################################################
save(distMat,beta, sigma2, phi, n,p,w, X_trainingdata, locations_trainingdata, Z_trainingdata, w_trainingdata, R_phi_trainingdata,X_testdata, locations_testdata, Z_testdata, w_testdata,
     file=paste0("G_phi",phiselect,"_Gen500_",sim,".RData")) 
  } 
}

################################################################################
################################################################################
