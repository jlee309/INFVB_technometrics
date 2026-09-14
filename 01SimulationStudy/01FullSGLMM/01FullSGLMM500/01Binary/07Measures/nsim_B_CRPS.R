###########################################################################################
###########################################################################################
# 11122024 
# CRPS Measure for the binary_nsim
rm(list=ls())
library(scoringRules)
getwd()
# setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/01Binary/07Measures")

################################################################################
################################################################################
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  
  crpsTotal_MCMC_mean_vector<-c()
  crpsTotal_INLA_mean_vector<-c()
  crpsTotal_VBphi_mean_vector<-c()
  crpsTotal_VBssqphi_mean_vector<-c()
  
  logsTotal_MCMC_median_vector<-c()
  logsTotal_INLA_meadian_vector<-c()
  logsTotal_VBphi_median_vector<-c()
  logsTotal_VBssqphi_median_vector<-c()
  
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name1 <- paste0("../01DataGeneration/B_phi", phiselect, "_Gen500_", sim, ".RData")
    file_name2 <- paste0("../02MCMC/B_phi", phiselect, "_Gen500_", sim, "_MCMC.RData")
    file_name3 <- paste0("../06INLA/B_phi", phiselect, "_Gen500_", sim, "_INLA.RData")
    file_name4 <- paste0("../03INFVBfixphi_ParallelComputing/B_phi", phiselect, "_Gen500_", sim, "_VBphiParallel.RData")
    file_name5 <- paste0("../08INFVBfixphi_SA/SA_B_phi", phiselect, "_Gen500_", sim, "_Full_INFVBphi.RData")
    file_name6 <- paste0("../09INFVBfixssqphi_SA/SA_B_phi", phiselect, "_Gen500_", sim, "_Full_INFVBssqphi.RData")
  
    load(file_name1)
    load(file_name2)
    load(file_name3)
    load(file_name4)
    load(file_name5)
    load(file_name6)
    
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)


# load(file="../01DataGeneration/B_phi1_Gen500_1.RData")
# load(file="../02MCMC/B_phi1_Gen500_1_MCMC.RData")
# load(file="../06INLA/B_phi1_Gen500_1_INLA.RData")
# load(file="../03INFVBfixphi_ParallelComputing/B_phi1_Gen500_1_VBphiParallel.RData")
# load(file="../08INFVBfixphi_SA/SA_B_phi1_Gen500_1_Full_INFVBphi.RData")
# load(file="../09INFVBfixssqphi_SA/SA_B_phi1_Gen500_1_Full_INFVBssqphi.RData")

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

# # Tail-area data
# tailData<-which(eta_testdata_real>quantile(eta_testdata_real, probs = 0.99))
# mean(crpsTotal_MCMC[tailData])
# median(crpsTotal_MCMC[tailData])
# median(logsTotal_MCMC[tailData])
# mean(logsTotal_MCMC[tailData])
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
  sampling_INLA<-rnorm(1000, mean=field_est_mean[hj], sd=field_est_sd[hj])
  crpsTotal_INLA[hj]<-crps_sample(y = eta_testdata_real[hj], # Continuous Ranked Probability Score
                                  dat = sampling_INLA,
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
# INFVBphi_testdata
###########################################################################################
###########################################################################################
crpsTotal_VBphi<-vector("numeric")
logsTotal_VBphi<-vector("numeric")

pt<-proc.time()
for(hj in 1:length(eta_testdata_real)){
  if(hj%%10==0){print(hj)}
  sampling_VBphi<-rnorm(1000, mean=mu_sample_pred_VB[hj], sd=sqrt(cov_sample_eta_testdata_VB_all[hj,hj]))
  crpsTotal_VBphi[hj]<-crps_sample(y = eta_testdata_real[hj], # Continuous Ranked Probability Score
                                   dat = sampling_VBphi,
                                   method = "edf")
  logsTotal_VBphi[hj]<-logs_sample(y = eta_testdata_real[hj], dat = sampling_VBphi) # Log Score
}
ptFinal2<-proc.time()-pt

# All of the data
print(ptFinal2)
median(crpsTotal_VBphi) 
mean(crpsTotal_VBphi)
median(logsTotal_VBphi)

###########################################################################################
###########################################################################################
# INFVBssqphi_testdata
###########################################################################################
###########################################################################################
crpsTotal_VBssqphi<-vector("numeric")
logsTotal_VBssqphi<-vector("numeric")


pt<-proc.time()
for(hj in 1:length(eta_testdata_real)){
  if(hj%%10==0){print(hj)}
  sampling_VBssqphi<-rnorm(1000, mean=mu_sample_pred_VBssqphi[hj], sd=sqrt(cov_sample_eta_testdata_VBssqphi_all[hj,hj]))
  crpsTotal_VBssqphi[hj]<-crps_sample(y = eta_testdata_real[hj], # Continuous Ranked Probability Score
                                      dat = sampling_VBssqphi,
                                      method = "edf")
  logsTotal_VBssqphi[hj]<-logs_sample(y = eta_testdata_real[hj], dat = sampling_VBssqphi) # Log Score
}
ptFinal2<-proc.time()-pt

# All of the data
print(ptFinal2)
median(crpsTotal_VBssqphi) 
mean(crpsTotal_VBssqphi)
median(logsTotal_VBssqphi)


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
crpsTotal_VBphi_mean<-mean(crpsTotal_VBphi)
crpsTotal_VBssqphi_mean<-mean(crpsTotal_VBssqphi)

logsTotal_MCMC_median<-median(logsTotal_MCMC)
logsTotal_INLA_meadian<-median(logsTotal_INLA)
logsTotal_VBphi_median<-median(logsTotal_VBphi)
logsTotal_VBssqphi_median<-median(logsTotal_VBssqphi)
###########################################################################################
###########################################################################################
crpsTotal_MCMC_mean_vector<-c(crpsTotal_MCMC_mean_vector,crpsTotal_MCMC_mean)
crpsTotal_INLA_mean_vector<-c(crpsTotal_INLA_mean_vector,crpsTotal_INLA_mean)
crpsTotal_VBphi_mean_vector<-c(crpsTotal_VBphi_mean_vector, crpsTotal_VBphi_mean)
crpsTotal_VBssqphi_mean_vector<-c(crpsTotal_VBssqphi_mean_vector, crpsTotal_VBssqphi_mean)

logsTotal_MCMC_median_vector<-c(logsTotal_MCMC_median_vector,logsTotal_MCMC_median)
logsTotal_INLA_meadian_vector<-c(logsTotal_INLA_meadian_vector,logsTotal_INLA_meadian)
logsTotal_VBphi_median_vector<-c(logsTotal_VBphi_median_vector, logsTotal_VBphi_median)
logsTotal_VBssqphi_median_vector<-c(logsTotal_VBssqphi_median_vector,logsTotal_VBssqphi_median)

  }
  save(crpsTotal_MCMC_mean_vector, crpsTotal_INLA_mean_vector, 
       crpsTotal_VBphi_mean_vector, crpsTotal_VBssqphi_mean_vector, 
       logsTotal_MCMC_median_vector, logsTotal_INLA_meadian_vector, 
       logsTotal_VBphi_median_vector, logsTotal_VBssqphi_median_vector, 
       file = paste0("B_phi",phiselect,"_Gen500_CRPS.RData"))
}

# rm(list=ls())
# load("B_phi2_Gen500_CRPS.RData")
# mean(crpsTotal_MCMC_mean_vector)
