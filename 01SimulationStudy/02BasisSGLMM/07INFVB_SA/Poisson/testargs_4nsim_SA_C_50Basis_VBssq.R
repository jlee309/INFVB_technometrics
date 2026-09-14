################################################################################
################################################################################
################################################################################
################################################################################
# 11032024 
# Advanced Sampling Method for INFVB_50bases
rm(list=ls())
library(mvtnorm);library(Matrix);
library(cli);library(viridis);library(fields);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm)
library(invgamma);library(pROC);library(Matrix);library(emulator)
library(argparser);library(tidyverse);

################################################################################
################################################################################
#Distributing our work to different CPU
args <- arg_parser('jin_test') %>%
  add_argument('--task', type = 'integer', help = 'Task number.') %>%
  parse_args()

# phi_test<-c(1,2,3,4) #it corresponds to c(0.1, 0.3, 0.5, 0.7)
phi_test<-4
nsim_test<-1:50 #should change this to 50 later
settings<-expand.grid(phi_test=phi_test,nsim_test=nsim_test)
phiselect<-settings$phi_test[args$task]
sim<-settings$nsim_test[args$task]


# Number of simulations
# n_simulations <- 1 #Number of simulations per each scenario
# phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed

# for(phiselect in 1:length(phiset) ){
#   cat("Running phiselect", phiselect, "\n")
  # for (sim in 1:n_simulations) {
  #   cat("Running simulation", sim, "\n")
    file_name1 <- paste0("../../01DataGeneration/phi",phiselect,"_Gen25k_",sim,"EigenBasis.RData")
    file_name2 <-paste0("../../01DataGeneration/phi",phiselect,"_Gen25k_",sim,"SpatialData.RData")
    file_name3 <-paste0("../../03INFVB_ParallelComputing/Poisson/50Basis_C_phi",phiselect,"_Gen25k_",sim,"_VBssqParallel.RData")
    
    load(file_name1)
    load(file_name2)
    load(file_name3)
    
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    

################################################################################
################################################################################
M<-refEigen_trainingdata[,1:numofbasis]
M_CV<-refEigen_testdata[,1:numofbasis]
dim(M);dim(M_CV)

################################################################################
################################################################################
# Sampling Method for estimating Eta and predicting Eta 
# By using Multivariate Normal for v(Beta, delta) 

################################################################################
# (1) Getting samples for Gamma(Beta, W)
################################################################################
numofsamples <- 50
sample_V_list <- lapply(1:iter.thetaphi, function(k) {
  if (k %% 100 == 0) print(k)
  thetaC_covMatgamma_Mat <- diag(thetaC_covMatgamma[k, ])
  rmvnorm(numofsamples, mean = thetaC_Meangamma[k, ], sigma = thetaC_covMatgamma_Mat)
})

################################################################################
# (2) Getting numofsamples Eta for each discretization
################################################################################
sample_eta_trainingdata_VB <- lapply(1:iter.thetaphi, function(i) {
  if (i %% 100 == 0) print(i)
  lapply(1:numofsamples, function(j) {
    sample_V <- sample_V_list[[i]][j, ]
    X_trainingdata %*% sample_V[1:p] + M %*% sample_V[(p + 1):(numofbasis + p)]
  })
})

sample_eta_testdata_VB <- lapply(1:iter.thetaphi, function(i) {
  if (i %% 100 == 0) print(i)
  lapply(1:numofsamples, function(j) {
    sample_V <- sample_V_list[[i]][j, ]
    X_testdata %*% sample_V[1:p] + M_CV %*% sample_V[(p + 1):(numofbasis + p)]
  })
})

################################################################################
# (3) Get the ELBO*Samples_Estimation
################################################################################
ELBOmean_list_est <- lapply(1:iter.thetaphi, function(i) {
  if (i %% 10 == 0) print(i)
  sapply(1:(0.8 * n), function(k) {
    mean_values <- sapply(1:numofsamples, function(j) {
      unlist(sample_eta_trainingdata_VB[[i]][j])[k]
    })
    mean(mean_values) * ELBOvectorweights[i]
  })
})

################################################################################
# (4) Making the sample covariance matrix for the estimation
################################################################################
cov_sample_eta_trainingdata_VB_all<-matrix(data=rep(0,(0.8*n)^2), nrow=(0.8*n), ncol=(0.8*n))
for(i in 1:iter.thetaphi){
  cov_sample_eta_trainingdata_VB<-c()
  for(j in 1:numofsamples){
    cov_sample_eta_trainingdata_VB<-
      cbind(cov_sample_eta_trainingdata_VB, unlist(sample_eta_trainingdata_VB[[i]][j]))
  }
  cov_sample_eta_trainingdata_VB<-t(cov_sample_eta_trainingdata_VB)
  cov_sample_eta_trainingdata_VB<-cov(cov_sample_eta_trainingdata_VB)
  assign(paste0("cov_sample_eta_trainingdata_VB_i",i),cov_sample_eta_trainingdata_VB*ELBOvectorweights[i] )
  dynamic_name<-paste0("cov_sample_eta_trainingdata_VB_i",i)
  cov_sample_eta_trainingdata_VB_i<-get(dynamic_name)
  cov_sample_eta_trainingdata_VB_all<-cov_sample_eta_trainingdata_VB_all+cov_sample_eta_trainingdata_VB_i
  # Optionally remove the dynamically created object after use
  rm(list = dynamic_name)
}

################################################################################
# (5) Prediction for ELBO*Samples_Prediction
################################################################################
ELBOmean_list_pred <- lapply(1:iter.thetaphi, function(i) {
  if (i %% 100 == 0) print(i)
  sapply(1:(0.2 * n), function(k) {
    mean_values <- sapply(1:numofsamples, function(j) {
      unlist(sample_eta_testdata_VB[[i]][j])[k]
    })
    mean(mean_values) * ELBOvectorweights[i]
  })
})

################################################################################
# (6) Making the sample covariance matrix for test data
################################################################################
cov_sample_eta_testdata_VB_all<-matrix(data=rep(0,(0.2*n)^2), nrow=(0.2*n), ncol=(0.2*n))
for(i in 1:iter.thetaphi){
  cov_sample_eta_testdata_VB<-c()
  for(j in 1:numofsamples){
    cov_sample_eta_testdata_VB<-
      cbind(cov_sample_eta_testdata_VB, unlist(sample_eta_testdata_VB[[i]][j]))
  }
  cov_sample_eta_testdata_VB<-t(cov_sample_eta_testdata_VB)
  cov_sample_eta_testdata_VB<-cov(cov_sample_eta_testdata_VB)
  assign(paste0("cov_sample_eta_testdata_VB_i",i),cov_sample_eta_testdata_VB*ELBOvectorweights[i] )
  dynamic_name<-paste0("cov_sample_eta_testdata_VB_i",i)
  cov_sample_eta_testdata_VB_i<-get(dynamic_name)
  cov_sample_eta_testdata_VB_all<-cov_sample_eta_testdata_VB_all+cov_sample_eta_testdata_VB_i
  # Optionally remove the dynamically created object after use
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-Delete for the memory#-#-#-#-#-#-#-#-#-#-#-#-#-#
  rm(list = dynamic_name)
}

################################################################################
# (7) Making mu_sample_est_VB
################################################################################
mu_sample_est_VB <- Reduce("+", ELBOmean_list_est)

################################################################################
# (8) Making mu_sample_pred_VB
################################################################################
mu_sample_pred_VB <- Reduce("+", ELBOmean_list_pred)

################################################################################
# (9) Clean up unnecessary objects from memory
################################################################################
rm(list = setdiff(ls(), c("cov_sample_eta_trainingdata_VB_all", 
                          "cov_sample_eta_testdata_VB_all", 
                          "mu_sample_est_VB", "mu_sample_pred_VB", 
                          "phiselect","sim", 
                          "n_simulations", "phiset", 
                          "set.seed","base_seed")))
################################################################################
# (10) Save only the required objects
################################################################################
save(cov_sample_eta_trainingdata_VB_all, cov_sample_eta_testdata_VB_all,
     mu_sample_est_VB, mu_sample_pred_VB,
     file = paste0("SA_50Basis_C_phi", phiselect, "_Gen25k_", sim, "_VBssqParallel.RData"))
