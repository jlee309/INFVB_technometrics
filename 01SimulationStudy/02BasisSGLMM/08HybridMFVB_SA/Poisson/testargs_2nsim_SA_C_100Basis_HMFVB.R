################################################################################
################################################################################
# 11032024 
# Count Data - Basis function MFVB Eta sampling method_100bases
################################################################################
################################################################################
rm(list=ls())
library(mvtnorm);library(Matrix);
library(cli);library(viridis);library(fields);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm)
library(invgamma);library(pROC);library(Matrix);library(emulator)
library(argparser);library(tidyverse);


#Distributing our work to different CPU
args <- arg_parser('jin_test') %>%
  add_argument('--task', type = 'integer', help = 'Task number.') %>%
  parse_args()

# phi_test<-c(1,2,3,4) #it corresponds to c(0.1, 0.3, 0.5, 0.7)
phi_test<-2
nsim_test<-1:50 #should change this to 50 later
settings<-expand.grid(phi_test=phi_test,nsim_test=nsim_test)
phiselect<-settings$phi_test[args$task]
sim<-settings$nsim_test[args$task]

set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


file_name1 <- paste0("../../01DataGeneration/phi",phiselect,"_Gen25k_",sim,"EigenBasis.RData")
file_name2 <-paste0("../../01DataGeneration/phi",phiselect,"_Gen25k_",sim,"SpatialData.RData")
file_name3 <-paste0("../../04HybridMFVB/Poisson/100Basis_C_phi",phiselect,"_Gen25k_",sim,"_HMFVB.RData")
file_name4 <-paste0("../../03INFVB_ParallelComputing/Poisson/100Basis_C_phi",phiselect,"_Gen25k_",sim,"_VBssqParallel.RData")

load(file_name1)
load(file_name2)
load(file_name3)
load(file_name4)

randomseed <- base_seed + sim + phiselect * 1000
set.seed(randomseed)


set.seed(12123)
################################################################################
################################################################################
M <- refEigen_trainingdata[, 1:numofbasis]
M_CV <- refEigen_testdata[, 1:numofbasis]

################################################################################
################################################################################
# Sampling Method for estimating Eta and predicting Eta 
# By using Multivariate Normal for v(Beta, W) 

################################################################################
# (1) Sampling Method for estimating Eta and predicting Eta
numofsamples <- 50
j <- nrow(MFVB_Meangamma)
MFVB_covMatgamma_Mat <- diag(MFVB_covMatgamma[j, ])
sample_V_MFVB <- rmvnorm(numofsamples, mean = MFVB_Meangamma[j, ], sigma = MFVB_covMatgamma_Mat)

################################################################################
# (2) Getting numofsamples Eta 
################################################################################
# (2) Getting numofsamples Eta
sample_eta_trainingdata_MFVB <- lapply(1:numofsamples, function(i) {
  X_trainingdata %*% sample_V_MFVB[i, c(1:p)] + M %*% sample_V_MFVB[i, c((p + 1):(numofbasis + p))]
})

sample_eta_testdata_MFVB <- lapply(1:numofsamples, function(i) {
  X_testdata %*% sample_V_MFVB[i, c(1:p)] + M_CV %*% sample_V_MFVB[i, c((p + 1):(numofbasis + p))]
})

################################################################################
# (3) Get the Mean of samples (note : Not a prediction now)
################################################################################
# (3) Get the Mean of samples
MFVB_mean_SA_Est <- sapply(1:(0.8 * n), function(k) {
  mean(sapply(1:numofsamples, function(i) sample_eta_trainingdata_MFVB[[i]][k]))
})



################################################################################
# (4) Making the sample covariance matrix for the estimation 
################################################################################
# (4) Making the sample covariance matrix for the estimation
cov_sample_eta_trainingdata_MFVB <- do.call(cbind, lapply(1:numofsamples, function(k) unlist(sample_eta_trainingdata_MFVB[[k]])))
cov_sample_eta_trainingdata_MFVB <- cov(t(cov_sample_eta_trainingdata_MFVB))

################################################################################
################################################################################
# Prediction 
################################################################################
# (3-1) Get Samples_Prediction
################################################################################
# Prediction
# (3-1) Get Samples_Prediction
MFVB_mean_SA_Pred <- sapply(1:(0.2 * n), function(k) {
  mean(sapply(1:numofsamples, function(i) sample_eta_testdata_MFVB[[i]][k]))
})


################################################################################
# (4-1) Making the sample covariance matrix for the estimation 
################################################################################
# (4-1) Making the sample covariance matrix for the estimation
cov_sample_eta_testdata_MFVB <- do.call(cbind, lapply(1:numofsamples, function(k) unlist(sample_eta_testdata_MFVB[[k]])))
cov_sample_eta_testdata_MFVB <- cov(t(cov_sample_eta_testdata_MFVB))


# ################################################################################
# # (6) Clean up unnecessary objects from memory
# ################################################################################
# rm(list = setdiff(ls(), c("MFVB_mean_SA_Est","cov_sample_eta_trainingdata_MFVB", 
#                           "MFVB_mean_SA_Pred","cov_sample_eta_testdata_MFVB")))

save(MFVB_mean_SA_Est,cov_sample_eta_trainingdata_MFVB, 
     MFVB_mean_SA_Pred,cov_sample_eta_testdata_MFVB, 
     file=paste0("SA_100Basis_C_phi",phiselect,"_Gen25k_",sim,"_HMFVB.RData")) 


