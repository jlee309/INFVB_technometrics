###########################################################################################
###########################################################################################
# 11162024 
# CRPS Measure for the binary_nsim_50Basis SGLMM
rm(list=ls())
library(scoringRules)
getwd()
# setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/02BasisSGLMM/11Measures/Poisson")
################################################################################
################################################################################
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


for(phiselect in 4:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  
  crpsTotal_MCMC_mean_vector<-c()
  crpsTotal_INLA_mean_vector<-c()
  crpsTotal_VBssq_mean_vector<-c()
  crpsTotal_HMFVB_mean_vector<-c()
  
  logsTotal_MCMC_median_vector<-c()
  logsTotal_INLA_meadian_vector<-c()
  logsTotal_VBssq_median_vector<-c()
  logsTotal_HMFVB_median_vector<-c()
  
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name1 <- paste0("../../01DataGeneration/phi", phiselect, "_Gen25k_", sim, "SpatialData.RData")
    file_name2 <- paste0("../../01DataGeneration/phi", phiselect, "_Gen25k_", sim, "EigenBasis.RData")
    file_name3 <- paste0("../../02MCMC/Poisson/50Basis_C_phi", phiselect, "_Gen25k_", sim, "_MCMC.RData")
    file_name4 <- paste0("../../03INFVB_ParallelComputing/Poisson/50Basis_C_phi", phiselect, "_Gen25k_", sim, "_VBssqParallel.RData")
    file_name5 <- paste0("../../07INFVB_SA/Poisson/SA_50Basis_C_phi", phiselect, "_Gen25k_", sim, "_VBssqParallel.RData")
    file_name6 <- paste0("../../08HybridMFVB_SA/Poisson/SA_50Basis_C_phi", phiselect, "_Gen25k_", sim, "_HMFVB.RData")    
    file_name7 <- paste0("../../06INLA/Poisson/50Basis_C_phi", phiselect, "_Gen25k_", sim, "_INLA.RData")    
    load(file_name1)
    load(file_name2)
    load(file_name3)
    load(file_name4)
    load(file_name5)
    load(file_name6)
    load(file_name7)
    
    
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    
    ###########################################################################################
    ###########################################################################################
    # MCMC_testdata
    ###########################################################################################
    ###########################################################################################
    eta_testdata_real <-X_testdata%*%beta+w_testdata
    crpsTotal_MCMC<-vector("numeric")
    logsTotal_MCMC<-vector("numeric")
    
    pt<-proc.time()
    for(hj in 1:length(eta_testdata_real)){
      if(hj%%10==0){print(hj)}
      crpsTotal_MCMC[hj]<-crps_sample(y = eta_testdata_real[hj], # Continuous Ranked Probability Score
                                      dat = eta_testdata_MCMC[hj,],
                                      method = "edf")
      logsTotal_MCMC[hj]<-logs_sample(y = eta_testdata_real[hj], dat = eta_testdata_MCMC[hj,]) # Log Score
    }
    ptFinal2<-proc.time()-pt
    
    
    # All of the data
    print(ptFinal2)
    median(crpsTotal_MCMC) 
    mean(crpsTotal_MCMC)
    median(logsTotal_MCMC)
    
    ###########################################################################################
    ###########################################################################################
    # INLA_testdata
    ###########################################################################################
    ###########################################################################################
    crpsTotal_INLA<-vector("numeric")
    logsTotal_INLA<-vector("numeric")
    
    
    pt<-proc.time()
    for(hj in 1:length(eta_testdata_real)){
      if(hj%%10==0){print(hj)}
      sampling_INLA<-rnorm(1000, mean=linear_predictor_mean[(0.8*n)+hj], sd=linear_predictor_sd[(0.8*n)+hj])
      crpsTotal_INLA[hj]<-crps_sample(y = eta_testdata_real[hj], # Continuous Ranked Probability Score
                                      dat = linear_predictor_mean,
                                      method = "edf")
      logsTotal_INLA[hj]<-logs_sample(y = eta_testdata_real[hj], dat = sampling_INLA) # Log Score
    }
    ptFinal2<-proc.time()-pt
    
    # All of the data
    print(ptFinal2)
    median(crpsTotal_INLA) 
    mean(crpsTotal_INLA)
    median(logsTotal_INLA)
    
    
    
    ###########################################################################################
    ###########################################################################################
    # INFVBssq_testdata
    ###########################################################################################
    ###########################################################################################
    crpsTotal_VBssq<-vector("numeric")
    logsTotal_VBssq<-vector("numeric")
    
    pt<-proc.time()
    for(hj in 1:length(eta_testdata_real)){
      if(hj%%10==0){print(hj)}
      sampling_VBssq<-rnorm(1000, mean=mu_sample_pred_VB[hj], sd=sqrt(cov_sample_eta_testdata_VB_all[hj,hj]))
      crpsTotal_VBssq[hj]<-crps_sample(y = eta_testdata_real[hj], # Continuous Ranked Probability Score
                                       dat = sampling_VBssq,
                                       method = "edf")
      logsTotal_VBssq[hj]<-logs_sample(y = eta_testdata_real[hj], dat = sampling_VBssq) # Log Score
    }
    ptFinal2<-proc.time()-pt
    
    # All of the data
    print(ptFinal2)
    median(crpsTotal_VBssq) 
    mean(crpsTotal_VBssq)
    median(logsTotal_VBssq)
    
    ###########################################################################################
    ###########################################################################################
    # HMFVB_testdata
    ###########################################################################################
    ###########################################################################################
    crpsTotal_HMFVB<-vector("numeric")
    logsTotal_HMFVB<-vector("numeric")
    
    
    pt<-proc.time()
    for(hj in 1:length(eta_testdata_real)){
      if(hj%%10==0){print(hj)}
      sampling_HMFVB<-rnorm(1000, mean=MFVB_mean_SA_Pred[hj], sd=sqrt(cov_sample_eta_testdata_MFVB[hj,hj]))
      crpsTotal_HMFVB[hj]<-crps_sample(y = eta_testdata_real[hj], # Continuous Ranked Probability Score
                                       dat = sampling_HMFVB,
                                       method = "edf")
      logsTotal_HMFVB[hj]<-logs_sample(y = eta_testdata_real[hj], dat = sampling_HMFVB) # Log Score
    }
    ptFinal2<-proc.time()-pt
    
    # All of the data
    print(ptFinal2)
    median(crpsTotal_HMFVB) 
    mean(crpsTotal_HMFVB)
    median(logsTotal_HMFVB)
    
    
    ###########################################################################################
    ###########################################################################################
    # # MCMC_Stan code 
    # # (note) Stan is only for 4 cases. (phi1_1, phi2_1, phi3_1, phi4_1) 
    # 
    # load(file="../10Stan/B_phi1_Gen500_1_MCMC_Stan.RData")
    # 
    # 
    # crpsTotal_MCMC_Stan<-vector("numeric")
    # logsTotal_MCMC_Stan<-vector("numeric")
    # 
    # pt<-proc.time()
    # for(hj in 1:length(eta_testdata_real)){
    #   if(hj%%10==0){print(hj)}
    #   crpsTotal_MCMC_Stan[hj]<-crps_sample(y = eta_testdata_real[hj], # Continuous Ranked Probability Score
    #                                        dat = eta_testdata_MCMC_Stan[hj,],
    #                                        method = "edf")
    #   logsTotal_MCMC_Stan[hj]<-logs_sample(y = eta_testdata_real[hj], dat = eta_testdata_MCMC_Stan[hj,]) # Log Score
    # }
    # ptFinal2<-proc.time()-pt
    # 
    # 
    # # All of the data
    # print(ptFinal2)
    # median(crpsTotal_MCMC_Stan) 
    # mean(crpsTotal_MCMC_Stan)
    # median(logsTotal_MCMC_Stan)
    
    ###########################################################################################
    ###########################################################################################
    crpsTotal_MCMC_mean<-mean(crpsTotal_MCMC)
    crpsTotal_INLA_mean<-mean(crpsTotal_INLA)
    crpsTotal_VBssq_mean<-mean(crpsTotal_VBssq)
    crpsTotal_HMFVB_mean<-mean(crpsTotal_HMFVB)
    
    logsTotal_MCMC_median<-median(logsTotal_MCMC)
    logsTotal_INLA_meadian<-median(logsTotal_INLA)
    logsTotal_VBssq_median<-median(logsTotal_VBssq)
    logsTotal_HMFVB_median<-median(logsTotal_HMFVB)
    ###########################################################################################
    ###########################################################################################
    crpsTotal_MCMC_mean_vector<-c(crpsTotal_MCMC_mean_vector,crpsTotal_MCMC_mean)
    crpsTotal_INLA_mean_vector<-c(crpsTotal_INLA_mean_vector,crpsTotal_INLA_mean)
    crpsTotal_VBssq_mean_vector<-c(crpsTotal_VBssq_mean_vector, crpsTotal_VBssq_mean)
    crpsTotal_HMFVB_mean_vector<-c(crpsTotal_HMFVB_mean_vector, crpsTotal_HMFVB_mean)
    
    logsTotal_MCMC_median_vector<-c(logsTotal_MCMC_median_vector,logsTotal_MCMC_median)
    logsTotal_INLA_meadian_vector<-c(logsTotal_INLA_meadian_vector,logsTotal_INLA_meadian)
    logsTotal_VBssq_median_vector<-c(logsTotal_VBssq_median_vector, logsTotal_VBssq_median)
    logsTotal_HMFVB_median_vector<-c(logsTotal_HMFVB_median_vector,logsTotal_HMFVB_median)
    
  }
  save(crpsTotal_MCMC_mean_vector, crpsTotal_INLA_mean_vector, 
       crpsTotal_VBssq_mean_vector, crpsTotal_HMFVB_mean_vector, 
       logsTotal_MCMC_median_vector, logsTotal_INLA_meadian_vector, 
       logsTotal_VBssq_median_vector, logsTotal_HMFVB_median_vector, 
       file = paste0("50Basis_C_phi",phiselect,"_Gen25k_CRPS.RData"))
}

