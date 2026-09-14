###########################################################################################
###########################################################################################
# 11122024 
# CRPS Measure for the Gaussian_nsim
rm(list=ls())
library(scoringRules)
getwd()
# setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/03Normal_Supplementary/07Measures")

################################################################################
################################################################################
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  
  crpsTotal_MCMC_mean_vector<-c()
  crpsTotal_VBphi_mean_vector<-c()
  
  logsTotal_MCMC_median_vector<-c()
  logsTotal_VBphi_median_vector<-c()
  
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


# load(file="../01DataGeneration/G_phi1_Gen500_1.RData")
# load(file="../02MCMC/G_phi1_Gen500_1_MCMC.RData")
# load(file="../06INLA/G_phi1_Gen500_1_INLA.RData")
# load(file="../03INFVBfixphi_ParallelComputing/G_phi1_Gen500_1_VBphiParallel.RData")
# load(file="../08INFVBfixphi_SA/SA_G_phi1_Gen500_1_Full_INFVBphi.RData")
# load(file="../09INFVBfixssqphi_SA/SA_G_phi1_Gen500_1_Full_INFVBssqphi.RData")

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
                                  dat = Z_prediction_MCMC[hj,],
                                  method = "edf")
  logsTotal_MCMC[hj]<-logs_sample(y = eta_testdata_real[hj], dat = Z_prediction_MCMC[hj,]) # Log Score
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
# # MCMC_Stan code 
# # (note) Stan is only for 4 cases. (phi1_1, phi2_1, phi3_1, phi4_1) 
# 
# load(file="../10Stan/G_phi1_Gen500_1_MCMC_Stan.RData")
# 
# 
# crpsTotal_MCMC_Stan<-vector("numeric")
# logsTotal_MCMC_Stan<-vector("numeric")
# 
# pt<-proc.time()
# for(hj in 1:length(eta_testdata_real)){
#   if(hj%%10==0){print(hj)}
#   crpsTotal_MCMC_Stan[hj]<-crps_sample(y = eta_testdata_real[hj], # Continuous Ranked Probability Score
#                                        dat = Z_prediction_MCMC_Stan[hj,],
#                                        method = "edf")
#   logsTotal_MCMC_Stan[hj]<-logs_sample(y = eta_testdata_real[hj], dat = Z_prediction_MCMC_Stan[hj,]) # Log Score
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
crpsTotal_VBphi_mean<-mean(crpsTotal_VBphi)

logsTotal_MCMC_median<-median(logsTotal_MCMC)
logsTotal_VBphi_median<-median(logsTotal_VBphi)
###########################################################################################
###########################################################################################
crpsTotal_MCMC_mean_vector<-c(crpsTotal_MCMC_mean_vector,crpsTotal_MCMC_mean)
crpsTotal_VBphi_mean_vector<-c(crpsTotal_VBphi_mean_vector, crpsTotal_VBphi_mean)

logsTotal_MCMC_median_vector<-c(logsTotal_MCMC_median_vector,logsTotal_MCMC_median)
logsTotal_VBphi_median_vector<-c(logsTotal_VBphi_median_vector, logsTotal_VBphi_median)

  }
  save(crpsTotal_MCMC_mean_vector, 
       crpsTotal_VBphi_mean_vector, 
       logsTotal_MCMC_median_vector, 
       logsTotal_VBphi_median_vector, 
       file = paste0("G_phi",phiselect,"_Gen500_CRPS.RData"))
}

# rm(list=ls())
# load("G_phi2_Gen500_CRPS.RData")
# mean(crpsTotal_MCMC_mean_vector)
