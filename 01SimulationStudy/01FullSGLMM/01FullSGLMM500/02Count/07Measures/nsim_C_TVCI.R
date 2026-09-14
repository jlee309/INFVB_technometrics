################################################################################
################################################################################
# True Value in Credible Interval _ Count data
################################################################################
################################################################################
# 11142024 
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
  
  TVCItrain_MCMC_mean_vector<-c()
  TVCItrain_INLA_mean_vector<-c()
  TVCItrain_VBphi_mean_vector<-c()
  TVCItrain_VBssqphi_mean_vector<-c()
  
  TVCItest_MCMC_mean_vector<-c()
  TVCItest_INLA_mean_vector<-c()
  TVCItest_VBphi_mean_vector<-c()
  TVCItest_VBssqphi_mean_vector<-c()
  
  
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name1 <- paste0("../01DataGeneration/C_phi", phiselect, "_Gen500_", sim, ".RData")
    file_name2 <- paste0("../02MCMC/C_phi", phiselect, "_Gen500_", sim, "_MCMC.RData")
    file_name3 <- paste0("../06INLA/C_phi", phiselect, "_Gen500_", sim, "_INLA.RData")
    file_name4 <- paste0("../03INFVBfixphi_ParallelComputing/C_phi", phiselect, "_Gen500_", sim, "_VBphiParallel.RData")
    file_name5 <- paste0("../08INFVBfixphi_SA/SA_C_phi", phiselect, "_Gen500_", sim, "_Full_INFVBphi.RData")
    file_name6 <- paste0("../09INFVBfixssqphi_SA/SA_C_phi", phiselect, "_Gen500_", sim, "_Full_INFVBssqphi.RData")
    
    load(file_name1)
    load(file_name2)
    load(file_name3)
    load(file_name4)
    load(file_name5)
    load(file_name6)
    
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
      lower_bound <- qnorm(0.025,mean = field_est_mean[k] , sd = field_est_sd[k])
      upper_bound <- qnorm(0.975,mean = field_est_mean[k] , sd = field_est_sd[k])
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
      lower_bound <- qnorm(0.025,mean = field_pred_mean[k] , sd = field_pred_sd[k])
      upper_bound <- qnorm(0.975,mean = field_pred_mean[k] , sd = field_pred_sd[k])
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
    TVCI_INFVBphi_train<-c()
    for(k in 1:(0.8*n)){
      # Assuming posterior_samples is your vector of posterior draws
      lower_bound <- qnorm(0.025,mean = mu_sample_est_VB[k] , sd = sqrt(cov_sample_eta_trainingdata_VB_all[k,k]))
      upper_bound <- qnorm(0.975,mean = mu_sample_est_VB[k] , sd = sqrt(cov_sample_eta_trainingdata_VB_all[k,k]))
      credible_interval <- c(lower_bound, upper_bound)
      TVCI_INFVBphi_train.cal<- (eta_trainingdata_real[k] >= credible_interval[1] && 
                                   eta_trainingdata_real[k] <= credible_interval[2])
      TVCI_INFVBphi_train<-c(TVCI_INFVBphi_train,TVCI_INFVBphi_train.cal)
    }
    sum(TVCI_INFVBphi_train)/length(TVCI_INFVBphi_train)
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
    TVCI_INFVBssqphi_train<-c()
    for(k in 1:(0.8*n)){
      # Assuming posterior_samples is your vector of posterior draws
      lower_bound <- qnorm(0.025,mean = mu_sample_est_VBssqphi[k] , sd = sqrt(cov_sample_eta_trainingdata_VBssqphi_all[k,k]))
      upper_bound <- qnorm(0.975,mean = mu_sample_est_VBssqphi[k] , sd = sqrt(cov_sample_eta_trainingdata_VBssqphi_all[k,k]))
      credible_interval <- c(lower_bound, upper_bound)
      TVCI_INFVBssqphi_train.cal<- (eta_trainingdata_real[k] >= credible_interval[1] && 
                                      eta_trainingdata_real[k] <= credible_interval[2])
      TVCI_INFVBssqphi_train<-c(TVCI_INFVBssqphi_train,TVCI_INFVBssqphi_train.cal)
    }
    sum(TVCI_INFVBssqphi_train)/length(TVCI_INFVBssqphi_train)
    ################################################################################
    ################################################################################
    #For test Data
    #For test Data
    TVCI_INFVBssqphi_test<-c()
    for(k in 1:(0.2*n)){
      # Assuming posterior_samples is your vector of posterior draws
      lower_bound <- qnorm(0.025,mean = mu_sample_pred_VBssqphi[k] , sd = sqrt(cov_sample_eta_testdata_VBssqphi_all[k,k]))
      upper_bound <- qnorm(0.975,mean = mu_sample_pred_VBssqphi[k] , sd = sqrt(cov_sample_eta_testdata_VBssqphi_all[k,k]))
      credible_interval <- c(lower_bound, upper_bound)
      TVCI_INFVBssqphi_test.cal<- (eta_testdata_real[k] >= credible_interval[1] && 
                                     eta_testdata_real[k] <= credible_interval[2])
      TVCI_INFVBssqphi_test<-c(TVCI_INFVBssqphi_test,TVCI_INFVBssqphi_test.cal)
    }
    sum(TVCI_INFVBssqphi_test)/length(TVCI_INFVBssqphi_test)
    
    
    
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    # Checking for eta_5
    ################################################################################
    # par(mfrow=c(1,1))
    # #MCMC
    # plot(density(eta_testdata_MCMC[5,]), main = "Density for Linear Predictor Prediction of test_Eta_5", col="red", ylim=c(0,1))
    # abline(v=eta_testdata_real[5,], col="green")
    # # INLA
    # xSeq<-seq(-3,3,length.out=iter.thetaphi)
    # ySeq<-dnorm(x = xSeq , mean = field_pred_mean[5] , sd = field_pred_sd[5])
    # lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)
    # 
    # 
    # # Sampling INFVBphi 
    # ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[5], sd=sqrt(cov_sample_eta_testdata_VB_all[5,5]))
    # lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="orange", lty=1)
    # abline(v=eta_testdata_real[5,], col="green")
    # lower_bound <- qnorm(0.025,mean = mu_sample_pred_VB[5] , sd = sqrt(cov_sample_eta_testdata_VB_all[5,5]))
    # upper_bound <- qnorm(0.975,mean = mu_sample_pred_VB[5] , sd = sqrt(cov_sample_eta_testdata_VB_all[5,5]))
    # abline(v=lower_bound, col="red")
    # abline(v=upper_bound, col="red")
    # 
    # 
    # 
    # # Sampling INFVBssqphi 
    # ySeq_VBfixssqphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VBssqphi[5], sd=sqrt(cov_sample_eta_testdata_VBssqphi_all[5,5]))
    # lines(x=xSeq , y = ySeq_VBfixssqphi_sampling , typ="l" , col="blue", lty=1)
    # 
    # # lower_bound <- qnorm(0.025,mean = mu_sample_pred_VBssqphi[5] , sd = sqrt(cov_sample_eta_testdata_VBssqphi_all[5,5]))
    # # upper_bound <- qnorm(0.975,mean = mu_sample_pred_VBssqphi[5] , sd = sqrt(cov_sample_eta_testdata_VBssqphi_all[5,5]))
    # # abline(v=lower_bound, col="red")
    # # abline(v=upper_bound, col="red")
    # 
    # 
    # legend("topright", lty=c(1,1,1), legend = c("INFVB1","INFVB2","MCMC", "INLA"), 
    #        col=c("orange","blue","red", "purple"))
    # 
    # 
    # 
    # ################################################################################
    # ################################################################################
    # # Checking for eta_16
    # ################################################################################
    # par(mfrow=c(1,1))
    # #MCMC
    # plot(density(eta_testdata_MCMC[16,]), main = "Density for Linear Predictor Prediction of test_Eta_16", col="red", ylim=c(0,1))
    # abline(v=eta_testdata_real[16,], col="green")
    # # INLA
    # xSeq<-seq(-3,3,length.out=iter.thetaphi)
    # ySeq<-dnorm(x = xSeq , mean = field_pred_mean[16] , sd = field_pred_sd[16])
    # lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)
    # 
    # 
    # # Sampling INFVBphi 
    # ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[16], sd=sqrt(cov_sample_eta_testdata_VB_all[16,16]))
    # lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="orange", lty=1)
    # abline(v=eta_testdata_real[16,], col="green")
    # # lower_bound <- qnorm(0.025,mean = mu_sample_pred_VB[16] , sd = sqrt(cov_sample_eta_testdata_VB_all[16,16]))
    # # upper_bound <- qnorm(0.975,mean = mu_sample_pred_VB[16] , sd = sqrt(cov_sample_eta_testdata_VB_all[16,16]))
    # # abline(v=lower_bound, col="purple")
    # # abline(v=upper_bound, col="purple")
    # 
    # 
    # 
    # # Sampling INFVBssqphi 
    # ySeq_VBfixssqphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VBssqphi[16], sd=sqrt(cov_sample_eta_testdata_VBssqphi_all[16,16]))
    # lines(x=xSeq , y = ySeq_VBfixssqphi_sampling , typ="l" , col="blue", lty=1)
    # 
    # # lower_bound <- qnorm(0.025,mean = mu_sample_pred_VBssqphi[16] , sd = sqrt(cov_sample_eta_testdata_VBssqphi_all[16,16]))
    # # upper_bound <- qnorm(0.975,mean = mu_sample_pred_VBssqphi[16] , sd = sqrt(cov_sample_eta_testdata_VBssqphi_all[16,16]))
    # # abline(v=lower_bound, col="red")
    # # abline(v=upper_bound, col="red")
    # 
    # 
    # legend("topright", lty=c(1,1,1), legend = c("INFVB1","INFVB2","MCMC", "INLA"), 
    #        col=c("orange","blue","red", "purple"))
    # 
    # 
    # 
    # ################################################################################
    # ################################################################################
    # # Checking for eta_19
    # ################################################################################
    # par(mfrow=c(1,1))
    # #MCMC
    # plot(density(eta_testdata_MCMC[19,]), main = "Density for Linear Predictor Prediction of test_Eta_19", col="red", ylim=c(0,1))
    # abline(v=eta_testdata_real[19,], col="green")
    # # INLA
    # xSeq<-seq(-3,3,length.out=iter.thetaphi)
    # ySeq<-dnorm(x = xSeq , mean = field_pred_mean[19] , sd = field_pred_sd[19])
    # lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)
    # 
    # 
    # # Sampling INFVBphi 
    # ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[19], sd=sqrt(cov_sample_eta_testdata_VB_all[19,19]))
    # lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="orange", lty=1)
    # abline(v=eta_testdata_real[19,], col="green")
    # lower_bound <- qnorm(0.025,mean = mu_sample_pred_VB[19] , sd = sqrt(cov_sample_eta_testdata_VB_all[19,19]))
    # upper_bound <- qnorm(0.975,mean = mu_sample_pred_VB[19] , sd = sqrt(cov_sample_eta_testdata_VB_all[19,19]))
    # abline(v=lower_bound, col="purple")
    # abline(v=upper_bound, col="purple")
    # 
    # 
    # 
    # # Sampling INFVBssqphi 
    # ySeq_VBfixssqphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VBssqphi[19], sd=sqrt(cov_sample_eta_testdata_VBssqphi_all[19,19]))
    # lines(x=xSeq , y = ySeq_VBfixssqphi_sampling , typ="l" , col="blue", lty=1)
    # 
    # lower_bound <- qnorm(0.025,mean = mu_sample_pred_VBssqphi[19] , sd = sqrt(cov_sample_eta_testdata_VBssqphi_all[19,19]))
    # upper_bound <- qnorm(0.975,mean = mu_sample_pred_VBssqphi[19] , sd = sqrt(cov_sample_eta_testdata_VBssqphi_all[19,19]))
    # abline(v=lower_bound, col="red")
    # abline(v=upper_bound, col="red")
    # 
    # 
    # legend("topright", lty=c(1,1,1), legend = c("INFVB1","INFVB2","MCMC", "INLA"), 
    #        col=c("orange","blue","red", "purple"))
    
    
    
    
    ###########################################################################################
    ###########################################################################################
    TVCItrain_MCMC_mean<-sum(TVCI_MCMC_train)/length(TVCI_MCMC_train)
    TVCItrain_INLA_mean<-sum(TVCI_INLA_train)/length(TVCI_INLA_train)
    TVCItrain_VBphi_mean<-sum(TVCI_INFVBphi_train)/length(TVCI_INFVBphi_train)
    TVCItrain_VBssqphi_mean<-sum(TVCI_INFVBssqphi_train)/length(TVCI_INFVBssqphi_train)
    
    TVCItest_MCMC_mean<-sum(TVCI_MCMC_test)/length(TVCI_MCMC_test)
    TVCItest_INLA_mean<-sum(TVCI_INLA_test)/length(TVCI_INLA_test)
    TVCItest_VBphi_mean<-sum(TVCI_INFVBphi_test)/length(TVCI_INFVBphi_test)
    TVCItest_VBssqphi_mean<-sum(TVCI_INFVBssqphi_test)/length(TVCI_INFVBssqphi_test)
    
    
    ###########################################################################################
    ###########################################################################################
    TVCItrain_MCMC_mean_vector<-c(TVCItrain_MCMC_mean_vector,TVCItrain_MCMC_mean)
    TVCItrain_INLA_mean_vector<-c(TVCItrain_INLA_mean_vector,TVCItrain_INLA_mean)
    TVCItrain_VBphi_mean_vector<-c(TVCItrain_VBphi_mean_vector, TVCItrain_VBphi_mean)
    TVCItrain_VBssqphi_mean_vector<-c(TVCItrain_VBssqphi_mean_vector, TVCItrain_VBssqphi_mean)
    
    TVCItest_MCMC_mean_vector<-c(TVCItest_MCMC_mean_vector,TVCItest_MCMC_mean)
    TVCItest_INLA_mean_vector<-c(TVCItest_INLA_mean_vector,TVCItest_INLA_mean)
    TVCItest_VBphi_mean_vector<-c(TVCItest_VBphi_mean_vector, TVCItest_VBphi_mean)
    TVCItest_VBssqphi_mean_vector<-c(TVCItest_VBssqphi_mean_vector, TVCItest_VBssqphi_mean)
    
  }
  save(TVCItrain_MCMC_mean_vector, TVCItrain_INLA_mean_vector, 
       TVCItrain_VBphi_mean_vector, TVCItrain_VBssqphi_mean_vector, 
       TVCItest_MCMC_mean_vector, TVCItest_INLA_mean_vector, 
       TVCItest_VBphi_mean_vector, TVCItest_VBssqphi_mean_vector, 
       file = paste0("C_phi",phiselect,"_Gen500_TVCI.RData"))
}

# getwd()
# 
# 
# 
# a1<-TVCItrain_MCMC_mean_vector
# a2<-TVCItrain_INLA_mean_vector 
# a3<-TVCItrain_VBphi_mean_vector 
# a4<-TVCItrain_VBssqphi_mean_vector 
# a5<-TVCItest_MCMC_mean_vector
# a6<-TVCItest_INLA_mean_vector 
# a7<-TVCItest_VBphi_mean_vector 
# a8<-TVCItest_VBssqphi_mean_vector 
# 
# TVCItrain_MCMC_mean_vector<-c(TVCItrain_MCMC_mean_vector,a1)
# TVCItrain_INLA_mean_vector<-c(TVCItrain_INLA_mean_vector,a2)
# TVCItrain_VBphi_mean_vector<-c(TVCItrain_VBphi_mean_vector,a3)
# TVCItrain_VBssqphi_mean_vector<-c(TVCItrain_VBssqphi_mean_vector,a4)
# 
# TVCItest_MCMC_mean_vector<-c(TVCItest_MCMC_mean_vector,a5)
# TVCItest_INLA_mean_vector<-c(TVCItest_INLA_mean_vector,a6)
# TVCItest_VBphi_mean_vector<-c(TVCItest_VBphi_mean_vector,a7)
# TVCItest_VBssqphi_mean_vector<-c(TVCItest_VBssqphi_mean_vector,a8)
# 
# save(TVCItrain_MCMC_mean_vector, TVCItrain_INLA_mean_vector, 
#      TVCItrain_VBphi_mean_vector, TVCItrain_VBssqphi_mean_vector, 
#      TVCItest_MCMC_mean_vector, TVCItest_INLA_mean_vector, 
#      TVCItest_VBphi_mean_vector, TVCItest_VBssqphi_mean_vector, 
#      file = "C_phi4_Gen500_TVCI.RData")



