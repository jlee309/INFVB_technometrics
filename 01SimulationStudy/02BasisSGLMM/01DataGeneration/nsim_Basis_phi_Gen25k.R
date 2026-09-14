################################################################################
################################################################################
# Data Generation for Binary and Poisson Data 
# using Basis functions n=25k observations 
# nsim
################################################################################
################################################################################
rm(list=ls())
# dev.off()
getwd()
library(fields) ; library(MASS) ; library(spam) ; library(RSpectra)
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(5153)
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    n=25000 # Sample Size

    
## Parameters
beta<-c(1,1) ; sigma2=1 ; phi=phiset[phiselect] ;
## Locations and Design Matrix (X)
X<-cbind(runif(n,-1,1), runif(n,-1,1))
locations<-cbind(runif(n),runif(n))
distMat<-as.matrix(rdist(locations)) # Distance Matrix
R_phi<-sigma2*exp(-distMat/phi) # Covariance Matrix using Exponential Covariance Function
################################################################################################
system.time(eigsUse<-eigs_sym(R_phi, k=300)) #k is the number of eigenvalues generates. 
set.seed(12345)
w<-eigsUse$vectors%*%(diag(sqrt(eigsUse$values))%*%rnorm(length(eigsUse$values))) 
XB<-X%*%beta
eta<-XB+w
lambda<-exp(eta)
summary(lambda)
prob<-exp(eta)/(1+exp(eta))
summary(prob)
obsPois<-rpois(n=n, lambda=lambda)
obsBin<-rbinom(n=n, size=1, prob=prob)


################################################################################
################################################################################
## Splitting into Two Data(Training, Test Data)
################################################################################
################################################################################
#Training Data (80% of the whole data)
modInd<-sample(1:n, n*0.8)
cvInd<-(1:n)[-modInd]
X_trainingdata<-X[modInd,]
locations_trainingdata<-locations[modInd,]
obsPois_trainingdata<-obsPois[modInd]
lambda_trainingdata<-lambda[modInd]
obsBin_trainingdata<-obsBin[modInd]
prob_trainingdata<-prob[modInd]
w_trainingdata<-w[modInd,]

#Test Data (20% of the whole data)
X_testdata<-X[cvInd,]
locations_testdata<-locations[cvInd,]
obsPois_testdata<-obsPois[cvInd]
lambda_testdata<-lambda[cvInd]
obsBin_testdata<-obsBin[cvInd]
prob_testdata<-prob[cvInd]
w_testdata<-w[cvInd,]


################################################################################
################################################################################
## Eigen Basis functions
################################################################################
################################################################################
## Basis Functions
## Making nX(number of basis) refEigen matrix
refEigen<-eigsUse$vectors
refEigen_trainingdata<-refEigen[modInd,]
refEigen_testdata<-refEigen[cvInd,]
################################################################################
################################################################################
save(refEigen_trainingdata,refEigen_testdata ,
      file=paste0("phi",phiselect,"_Gen25k_",sim,"EigenBasis.RData"))


save(X_trainingdata , locations_trainingdata , obsPois_trainingdata, lambda_trainingdata , 
     obsBin_trainingdata, prob_trainingdata , w_trainingdata , 
     X_testdata , locations_testdata , obsPois_testdata , lambda_testdata , 
     obsBin_testdata , prob_testdata , w_testdata , 
     beta , sigma2 , phi , n , prob,w,
     file=paste0("phi",phiselect,"_Gen25k_",sim,"SpatialData.RData"))
  } 
}