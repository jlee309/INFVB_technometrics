################################################################################
################################################################################
# True Value in Credible Interval 50basis
################################################################################
################################################################################
# 11142024 
################################################################################
################################################################################
getwd()
rm(list=ls())
# setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/02BasisSGLMM/11Measures/Binary")
################################################################################
################################################################################
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  
  TVCItrain_MCMC_mean_vector<-c()
  TVCItrain_INLA_mean_vector<-c()
  TVCItrain_VBssq_mean_vector<-c()
  TVCItrain_HMFVB_mean_vector<-c()
  
  TVCItest_MCMC_mean_vector<-c()
  TVCItest_INLA_mean_vector<-c()
  TVCItest_VBssq_mean_vector<-c()
  TVCItest_HMFVB_mean_vector<-c()
  
  
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name1 <- paste0("../../01DataGeneration/phi", phiselect, "_Gen25k_", sim, "SpatialData.RData")
    file_name2 <- paste0("../../01DataGeneration/phi", phiselect, "_Gen25k_", sim, "EigenBasis.RData")
    file_name3 <- paste0("../../02MCMC/Binary/50Basis_B_phi", phiselect, "_Gen25k_", sim, "_MCMC.RData")
    file_name4 <- paste0("../../03INFVB_ParallelComputing/Binary/50Basis_B_phi", phiselect, "_Gen25k_", sim, "_VBssqParallel.RData")
    file_name5 <- paste0("../../07INFVB_SA/Binary/SA_50Basis_B_phi", phiselect, "_Gen25k_", sim, "_VBssqParallel.RData")
    file_name6 <- paste0("../../08HybridMFVB_SA/Binary/SA_50Basis_B_phi", phiselect, "_Gen25k_", sim, "_HMFVB.RData")    
    file_name7 <- paste0("../../06INLA/Binary/50Basis_B_phi", phiselect, "_Gen25k_", sim, "_INLA.RData")    
    load(file_name1)
    load(file_name2)
    load(file_name3)
    load(file_name4)
    load(file_name5)
    load(file_name6)
    load(file_name7)
    
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    
    # The real Data
    eta_trainingdata_real <-X_trainingdata%*%beta+w_trainingdata
    eta_testdata_real <-X_testdata%*%beta+w_testdata
    
    
    ################################################################################
    ################################################################################
    # (1) MCMC methods 
    #For training Data
    TVCI_MCMC_train<-c()
    for(k in 1:(0.8*n)){
      posterior_samples<-eta_trainingdata_MCMC[k,]
      # Assuming posterior_samples is your vector of posterior draws
      lower_bound <- quantile(posterior_samples, 0.025)
      upper_bound <- quantile(posterior_samples, 0.975)
      credible_interval <- c(lower_bound, upper_bound)
      TVCI_MCMC_train.cal<- (eta_trainingdata_real[k] >= credible_interval[1] && 
                               eta_trainingdata_real[k] <= credible_interval[2])
      TVCI_MCMC_train<-c(TVCI_MCMC_train,TVCI_MCMC_train.cal)
    }
    sum(TVCI_MCMC_train)/length(TVCI_MCMC_train)
    
    ################################################################################
    ################################################################################
    #For test Data
    TVCI_MCMC_test<-c()
    for(k in 1:(0.2*n)){
      posterior_samples<-eta_testdata_MCMC[k,]
      # Assuming posterior_samples is your vector of posterior draws
      lower_bound <- quantile(posterior_samples, 0.025)
      upper_bound <- quantile(posterior_samples, 0.975)
      credible_interval <- c(lower_bound, upper_bound)
      TVCI_MCMC_test.cal<- (eta_testdata_real[k] >= credible_interval[1] && 
                              eta_testdata_real[k] <= credible_interval[2])
      TVCI_MCMC_test<-c(TVCI_MCMC_test,TVCI_MCMC_test.cal)
    }
    sum(TVCI_MCMC_test)/length(TVCI_MCMC_test)
    
    
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    # True Value in Credible Interval (2) INLA methods 
    ################################################################################
    ################################################################################
    #For training Data
    TVCI_INLA_train<-c()
    for(k in 1:(0.8*n)){
      # Assuming posterior_samples is your vector of posterior draws
      lower_bound <- qnorm(0.025,mean = linear_predictor_mean[k] , sd = linear_predictor_sd[k])
      upper_bound <- qnorm(0.975,mean = linear_predictor_mean[k] , sd = linear_predictor_sd[k])
      credible_interval <- c(lower_bound, upper_bound)
      TVCI_INLA_train.cal<- (eta_trainingdata_real[k] >= credible_interval[1] && 
                               eta_trainingdata_real[k] <= credible_interval[2])
      TVCI_INLA_train<-c(TVCI_INLA_train,TVCI_INLA_train.cal)
    }
    sum(TVCI_INLA_train)/length(TVCI_INLA_train)
    ################################################################################
    ################################################################################
    #For test Data
    #For test Data
    TVCI_INLA_test<-c()
    for(k in 1:(0.2*n)){
      # Assuming posterior_samples is your vector of posterior draws
      lower_bound <- qnorm(0.025,mean = linear_predictor_mean[(0.8*n)+k] , sd = linear_predictor_sd[(0.8*n)+k])
      upper_bound <- qnorm(0.975,mean = linear_predictor_mean[(0.8*n)+k] , sd = linear_predictor_sd[(0.8*n)+k])
      credible_interval <- c(lower_bound, upper_bound)
      TVCI_INLA_test.cal<- (eta_testdata_real[k] >= credible_interval[1] && 
                              eta_testdata_real[k] <= credible_interval[2])
      TVCI_INLA_test<-c(TVCI_INLA_test,TVCI_INLA_test.cal)
    }
    sum(TVCI_INLA_test)/length(TVCI_INLA_test)
    
    
    
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    # True Value in Credible Interval (3) INFVBfixphi method 
    ################################################################################
    ################################################################################
    #For training Data
    TVCI_INFVBssq_train<-c()
    for(k in 1:(0.8*n)){
      # Assuming posterior_samples is your vector of posterior draws
      lower_bound <- qnorm(0.025,mean = mu_sample_est_VB[k] , sd = sqrt(cov_sample_eta_trainingdata_VB_all[k,k]))
      upper_bound <- qnorm(0.975,mean = mu_sample_est_VB[k] , sd = sqrt(cov_sample_eta_trainingdata_VB_all[k,k]))
      credible_interval <- c(lower_bound, upper_bound)
      TVCI_INFVBssq_train.cal<- (eta_trainingdata_real[k] >= credible_interval[1] && 
                                   eta_trainingdata_real[k] <= credible_interval[2])
      TVCI_INFVBssq_train<-c(TVCI_INFVBssq_train,TVCI_INFVBssq_train.cal)
    }
    sum(TVCI_INFVBssq_train)/length(TVCI_INFVBssq_train)
    ################################################################################
    ################################################################################
    #For test Data
    TVCI_INFVBssq_test<-c()
    for(k in 1:(0.2*n)){
      # Assuming posterior_samples is your vector of posterior draws
      lower_bound <- qnorm(0.025,mean = mu_sample_pred_VB[k] , sd = sqrt(cov_sample_eta_testdata_VB_all[k,k]))
      upper_bound <- qnorm(0.975,mean = mu_sample_pred_VB[k] , sd = sqrt(cov_sample_eta_testdata_VB_all[k,k]))
      credible_interval <- c(lower_bound, upper_bound)
      TVCI_INFVBssq_test.cal<- (eta_testdata_real[k] >= credible_interval[1] && 
                                  eta_testdata_real[k] <= credible_interval[2])
      TVCI_INFVBssq_test<-c(TVCI_INFVBssq_test,TVCI_INFVBssq_test.cal)
    }
    sum(TVCI_INFVBssq_test)/length(TVCI_INFVBssq_test)
    
    
    
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    # True Value in Credible Interval (3) INFVBfixssqphi method 
    ################################################################################
    ################################################################################
    #For training Data
    TVCI_INFHMFVB_train<-c()
    for(k in 1:(0.8*n)){
      # Assuming posterior_samples is your vector of posterior draws
      lower_bound <- qnorm(0.025,mean = MFVB_mean_SA_Est[k] , sd = sqrt(cov_sample_eta_trainingdata_MFVB[k,k]))
      upper_bound <- qnorm(0.975,mean = MFVB_mean_SA_Est[k] , sd = sqrt(cov_sample_eta_trainingdata_MFVB[k,k]))
      credible_interval <- c(lower_bound, upper_bound)
      TVCI_INFHMFVB_train.cal<- (eta_trainingdata_real[k] >= credible_interval[1] && 
                                   eta_trainingdata_real[k] <= credible_interval[2])
      TVCI_INFHMFVB_train<-c(TVCI_INFHMFVB_train,TVCI_INFHMFVB_train.cal)
    }
    sum(TVCI_INFHMFVB_train)/length(TVCI_INFHMFVB_train)
    ################################################################################
    ################################################################################
    #For test Data
    #For test Data
    TVCI_INFHMFVB_test<-c()
    for(k in 1:(0.2*n)){
      # Assuming posterior_samples is your vector of posterior draws
      lower_bound <- qnorm(0.025,mean = MFVB_mean_SA_Pred[k] , sd = sqrt(cov_sample_eta_testdata_MFVB[k,k]))
      upper_bound <- qnorm(0.975,mean = MFVB_mean_SA_Pred[k] , sd = sqrt(cov_sample_eta_testdata_MFVB[k,k]))
      credible_interval <- c(lower_bound, upper_bound)
      TVCI_INFHMFVB_test.cal<- (eta_testdata_real[k] >= credible_interval[1] && 
                                  eta_testdata_real[k] <= credible_interval[2])
      TVCI_INFHMFVB_test<-c(TVCI_INFHMFVB_test,TVCI_INFHMFVB_test.cal)
    }
    sum(TVCI_INFHMFVB_test)/length(TVCI_INFHMFVB_test)
    
    
    
    
    
    ###########################################################################################
    ###########################################################################################
    TVCItrain_MCMC_mean<-sum(TVCI_MCMC_train)/length(TVCI_MCMC_train)
    TVCItrain_INLA_mean<-sum(TVCI_INLA_train)/length(TVCI_INLA_train)
    TVCItrain_VBssq_mean<-sum(TVCI_INFVBssq_train)/length(TVCI_INFVBssq_train)
    TVCItrain_HMFVB_mean<-sum(TVCI_INFHMFVB_train)/length(TVCI_INFHMFVB_train)
    
    TVCItest_MCMC_mean<-sum(TVCI_MCMC_test)/length(TVCI_MCMC_test)
    TVCItest_INLA_mean<-sum(TVCI_INLA_test)/length(TVCI_INLA_test)
    TVCItest_VBssq_mean<-sum(TVCI_INFVBssq_test)/length(TVCI_INFVBssq_test)
    TVCItest_HMFVB_mean<-sum(TVCI_INFHMFVB_test)/length(TVCI_INFHMFVB_test)
    
    
    ###########################################################################################
    ###########################################################################################
    TVCItrain_MCMC_mean_vector<-c(TVCItrain_MCMC_mean_vector,TVCItrain_MCMC_mean)
    TVCItrain_INLA_mean_vector<-c(TVCItrain_INLA_mean_vector,TVCItrain_INLA_mean)
    TVCItrain_VBssq_mean_vector<-c(TVCItrain_VBssq_mean_vector, TVCItrain_VBssq_mean)
    TVCItrain_HMFVB_mean_vector<-c(TVCItrain_HMFVB_mean_vector, TVCItrain_HMFVB_mean)
    
    TVCItest_MCMC_mean_vector<-c(TVCItest_MCMC_mean_vector,TVCItest_MCMC_mean)
    TVCItest_INLA_mean_vector<-c(TVCItest_INLA_mean_vector,TVCItest_INLA_mean)
    TVCItest_VBssq_mean_vector<-c(TVCItest_VBssq_mean_vector, TVCItest_VBssq_mean)
    TVCItest_HMFVB_mean_vector<-c(TVCItest_HMFVB_mean_vector, TVCItest_HMFVB_mean)
    
  }
  save(TVCItrain_MCMC_mean_vector, TVCItrain_INLA_mean_vector, 
       TVCItrain_VBssq_mean_vector, TVCItrain_HMFVB_mean_vector, 
       TVCItest_MCMC_mean_vector, TVCItest_INLA_mean_vector, 
       TVCItest_VBssq_mean_vector, TVCItest_HMFVB_mean_vector, 
       file = paste0("50Basis_B_phi",phiselect,"_Gen25k_TVCI.RData"))
}

