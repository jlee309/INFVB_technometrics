###########################################################################################
###########################################################################################
# 11122024 
# CRPS Measure for the binary 
###########################################################################################
###########################################################################################
# MCMC_testdata
###########################################################################################
###########################################################################################
rm(list=ls())
getwd()
setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/01Binary/07Measures")
library(scoringRules)
load(file="../01DataGeneration/B_phi1_Gen500_1.RData")
load(file="../02MCMC/B_phi1_Gen500_1_MCMC.RData")

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
load(file="../06INLA/B_phi1_Gen500_1_INLA.RData")

crpsTotal_INLA<-vector("numeric")
logsTotal_INLA<-vector("numeric")

set.seed(124)
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
load(file="../03INFVBfixphi_ParallelComputing/B_phi1_Gen500_1_VBphiParallel.RData")
load(file="../08INFVBfixphi_SA/SA_B_phi1_Gen500_1_Full_INFVBphi.RData")

crpsTotal_VBphi<-vector("numeric")
logsTotal_VBphi<-vector("numeric")

set.seed(124)
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
load(file="../09INFVBfixssqphi_SA/SA_B_phi1_Gen500_1_Full_INFVBssqphi.RData")
crpsTotal_VBssqphi<-vector("numeric")
logsTotal_VBssqphi<-vector("numeric")

set.seed(124)
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
# MCMC_Stan code 
# (note) Stan is only for 4 cases. (phi1_1, phi2_1, phi3_1, phi4_1) 
load(file="../10Stan/B_phi1_Gen500_1_MCMC_Stan.RData")


crpsTotal_MCMC_Stan<-vector("numeric")
logsTotal_MCMC_Stan<-vector("numeric")

pt<-proc.time()
for(hj in 1:length(eta_testdata_real)){
  if(hj%%10==0){print(hj)}
  crpsTotal_MCMC_Stan[hj]<-crps_sample(y = eta_testdata_real[hj], # Continuous Ranked Probability Score
                                  dat = eta_testdata_MCMC_Stan[hj,],
                                  method = "edf")
  logsTotal_MCMC_Stan[hj]<-logs_sample(y = eta_testdata_real[hj], dat = eta_testdata_MCMC_Stan[hj,]) # Log Score
}
ptFinal2<-proc.time()-pt


# All of the data
print(ptFinal2)
median(crpsTotal_MCMC_Stan) 
mean(crpsTotal_MCMC_Stan)
median(logsTotal_MCMC_Stan)

###########################################################################################
###########################################################################################
# Checking
#
#
mean(crpsTotal_MCMC)
mean(crpsTotal_INLA)
mean(crpsTotal_VBphi)
mean(crpsTotal_VBssqphi)
mean(crpsTotal_MCMC_Stan)

median(logsTotal_MCMC)
median(logsTotal_INLA)
median(logsTotal_VBphi)
median(logsTotal_VBssqphi)
median(logsTotal_MCMC_Stan)
###########################################################################################
###########################################################################################


aucVal
aucVal_INLA
aucVal_Stan
aucVal_VB

