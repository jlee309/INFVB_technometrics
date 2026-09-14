################################################################################
################################################################################
# True Value in Credible Interval 
################################################################################
################################################################################
# 11162024 
################################################################################
################################################################################
getwd()
rm(list=ls())
################################################################################
################################################################################
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  
  TVCItest_MCMC_mean_vector<-c()
  TVCItest_VBphi_mean_vector<-c()
  
  
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name1 <- paste0("../01DataGeneration/G_phi", phiselect, "_Gen500_", sim, ".RData")
    file_name2 <- paste0("../02MCMC/G_phi", phiselect, "_Gen500_", sim, "_MCMC.RData")
    file_name3 <- paste0("../03INFVBfixphi_ParallelComputing/G_phi", phiselect, "_Gen500_", sim, "_VBphiParallel.RData")
    file_name4 <- paste0("../05INFVBfixphi_SA/SA_G_phi", phiselect, "_Gen500_", sim, "_Full_INFVBphi.RData")
    
    load(file_name1)
    load(file_name2)
    load(file_name3)
    load(file_name4)
    
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    
    # The real Data
    # eta_trainingdata_real <-X_trainingdata%*%beta+w_trainingdata
    # eta_testdata_real <-X_testdata%*%beta+w_testdata

    Sigma_11 <- sigma2 * exp(-distMat[c(1:(0.8 * n)), c(1:(0.8 * n))] / sigma2)
    chol_inv_Sigma_11 <- chol(Sigma_11)
    inv_Sigma_11 <- chol2inv(chol_inv_Sigma_11)
    Sigma_21 <- sigma2* exp(-distMat[c((0.8 * n + 1):n), c(1:(0.8 * n))] / sigma2)
    eta_testdata_real<-X_testdata%*%beta+Sigma_21%*%inv_Sigma_11%*%(Z_trainingdata-X_trainingdata%*%beta)
    
    
    ################################################################################
    ################################################################################
    # (1) MCMC methods 
    ################################################################################
    ################################################################################
    #For test Data
    TVCI_MCMC_test<-c()
    for(k in 1:(0.2*n)){
      posterior_samples<-Z_prediction_MCMC[k,]
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
    # True Value in Credible Interval (2) INFVBfixphi method 
    ################################################################################
    ################################################################################
    #For test Data
    TVCI_INFVBphi_test<-c()
    for(k in 1:(0.2*n)){
      # Assuming posterior_samples is your vector of posterior draws
      lower_bound <- qnorm(0.025,mean = mu_sample_pred_VB[k] , sd = sqrt(cov_sample_eta_testdata_VB_all[k,k]))
      upper_bound <- qnorm(0.975,mean = mu_sample_pred_VB[k] , sd = sqrt(cov_sample_eta_testdata_VB_all[k,k]))
      credible_interval <- c(lower_bound, upper_bound)
      TVCI_INFVBphi_test.cal<- (eta_testdata_real[k] >= credible_interval[1] && 
                                  eta_testdata_real[k] <= credible_interval[2])
      TVCI_INFVBphi_test<-c(TVCI_INFVBphi_test,TVCI_INFVBphi_test.cal)
    }
    sum(TVCI_INFVBphi_test)/length(TVCI_INFVBphi_test)
    
    
    
    
    ###########################################################################################
    ###########################################################################################
    TVCItest_MCMC_mean<-sum(TVCI_MCMC_test)/length(TVCI_MCMC_test)
    TVCItest_VBphi_mean<-sum(TVCI_INFVBphi_test)/length(TVCI_INFVBphi_test)
    
    
    ###########################################################################################
    ###########################################################################################
    TVCItest_MCMC_mean_vector<-c(TVCItest_MCMC_mean_vector,TVCItest_MCMC_mean)
    TVCItest_VBphi_mean_vector<-c(TVCItest_VBphi_mean_vector, TVCItest_VBphi_mean)
  }
  save(TVCItest_MCMC_mean_vector, 
       TVCItest_VBphi_mean_vector, 
       file = paste0("G_phi",phiselect,"_Gen500_TVCI.RData"))
}
