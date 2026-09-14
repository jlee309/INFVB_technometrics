setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/01Binary/05Summary_Figures")
getwd()

rm(list=ls())
library(mvtnorm)
library(cli);library(viridis);library(fields);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm)
library(invgamma);library(pROC);library(Matrix);library(emulator)


load(file="../01DataGeneration/B_phi1_Gen500_1.RData")
load(file="../02MCMC/B_phi1_Gen500_1_MCMC.RData")
load(file="../03INFVBfixphi_ParallelComputing/B_phi1_Gen500_1_VBphiParallel.RData")
load(file="../08INFVBfixphi_SA/SA_B_phi1_Gen500_1_Full_INFVBphi.RData")
load(file="../09INFVBfixssqphi_SA/SA_B_phi1_Gen500_1_Full_INFVBssqphi.RData")
load(file="../06INLA/B_phi1_Gen500_1_INLA.RData")
load(file="../10Stan/B_phi1_Gen500_1_MCMC_Stan.RData")

dev.off()
################################################################################
# (5) Plotting Etas in the training datasets
################################################################################
# The real Data
eta_trainingdata_real <-X_trainingdata%*%beta+w_trainingdata
eta_testdata_real <-X_testdata%*%beta+w_testdata
################################################################################
# Plot one of the Eta_3_estimation
################################################################################
par(mfrow=c(1,1))
# Eta_3
#MCMC
plot(density(eta_trainingdata_MCMC[3,]), main = "Density for Linear Predictor Eta_3", col="red", ylim=c(0,1.5))
abline(v=eta_trainingdata_real[3,], col="green")
# # INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_est_mean[3] , sd = field_est_sd[3])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

# Sampling INFVBphi
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VB[3], sd=sqrt(cov_sample_eta_trainingdata_VB_all[3,3]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="orange", lty=1)

# Sampling INFVBssqphi
ySeq_VBfixssqphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VBssqphi[3], sd=sqrt(cov_sample_eta_trainingdata_VBssqphi_all[3,3]))
lines(x=xSeq , y = ySeq_VBfixssqphi_sampling , typ="l" , col="blue", lty=1)
# Stan
lines(density(eta_trainingdata_MCMC_Stan[3,]))


legend("topright", lty=c(1,1,1), legend = c("INFVB1","INFVB2","MCMC", "INLA","STAN"), 
       col=c("orange","blue","red", "purple","black"))


################################################################################
# Plot one of the Eta_10_estimation
################################################################################
# Eta_10
#MCMC
plot(density(eta_trainingdata_MCMC[10,]), main = "Density for Linear Predictor Eta_10", col="red", ylim=c(0,2))
abline(v=eta_trainingdata_real[10,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_est_mean[10] , sd = field_est_sd[10])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

# Sampling INFVBphi
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VB[10], sd=sqrt(cov_sample_eta_trainingdata_VB_all[10,10]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="orange", lty=1)

# Sampling INFVBssqphi
ySeq_VBfixssqphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VBssqphi[10], sd=sqrt(cov_sample_eta_trainingdata_VBssqphi_all[10,10]))
lines(x=xSeq , y = ySeq_VBfixssqphi_sampling , typ="l" , col="blue", lty=1)

# Stan
lines(density(eta_trainingdata_MCMC_Stan[10,]))

legend("topright", lty=c(1,1,1), legend = c("INFVB1","INFVB2","MCMC", "INLA","STAN"), 
       col=c("orange","blue","red", "purple","black"))


################################################################################
# Plot one of the Eta_25_estimation
################################################################################
# Eta_25
#MCMC
plot(density(eta_trainingdata_MCMC[25,]), main = "Density for Linear Predictor Eta_25", col="red", ylim=c(0,2))
abline(v=eta_trainingdata_real[25,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_est_mean[25] , sd = field_est_sd[25])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

# Sampling INFVBphi
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VB[25], sd=sqrt(cov_sample_eta_trainingdata_VB_all[25,25]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="orange", lty=1)

# Sampling INFVBssqphi
ySeq_VBfixssqphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VBssqphi[25], sd=sqrt(cov_sample_eta_trainingdata_VBssqphi_all[25,25]))
lines(x=xSeq , y = ySeq_VBfixssqphi_sampling , typ="l" , col="blue", lty=1)

# INLA
lines(density(eta_trainingdata_MCMC_Stan[25,]))

legend("topright", lty=c(1,1,1), legend = c("INFVB1","INFVB2","MCMC", "INLA","STAN"), 
       col=c("orange","blue","red", "purple","black"))

# legend("topright", lty=c(1,1,1), legend = c("MCMC", "INLA"), 
#        col=c("red", "purple"))

################################################################################
# Plot one of the Eta_40_estimation
################################################################################
# Eta_40
#MCMC
plot(density(eta_trainingdata_MCMC[40,]), main = "Density for Linear Predictor Eta_40", col="red", ylim=c(0,1.5))
abline(v=eta_trainingdata_real[40,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_est_mean[40] , sd = field_est_sd[40])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

# # # # # test for sampling method for eta_40
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VB[40], sd=sqrt(cov_sample_eta_trainingdata_VB_all[40,40]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="orange", lty=1)

# Sampling INFVBssqphi
ySeq_VBfixssqphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VBssqphi[40], sd=sqrt(cov_sample_eta_trainingdata_VBssqphi_all[40,40]))
lines(x=xSeq , y = ySeq_VBfixssqphi_sampling , typ="l" , col="blue", lty=1)


# INLA
lines(density(eta_trainingdata_MCMC_Stan[40,]))

legend("topright", lty=c(1,1,1), legend = c("INFVB1","INFVB2","MCMC", "INLA","STAN"), 
       col=c("orange","blue","red", "purple","black"))


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
################################################################################
# (5-1) Plotting Etas in the test datasets
################################################################################
################################################################################
# Plot one of the Eta_3_Prediction
################################################################################
#MCMC
plot(density(eta_testdata_MCMC[3,]), main = "Density for Linear Predictor Prediction of test_Eta_3", col="red", ylim=c(0,2))
abline(v=eta_testdata_real[3,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_pred_mean[3] , sd = field_pred_sd[3])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

# Sampling INFVBphi 
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[3], sd=sqrt(cov_sample_eta_testdata_VB_all[3,3]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="orange", lty=1)

# Sampling INFVBssqphi 
ySeq_VBfixssqphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VBssqphi[3], sd=sqrt(cov_sample_eta_testdata_VBssqphi_all[3,3]))
lines(x=xSeq , y = ySeq_VBfixssqphi_sampling , typ="l" , col="blue", lty=1)


# Stan
lines(density(eta_testdata_MCMC_Stan[3,]))

legend("topright", lty=c(1,1,1), legend = c("INFVB1","INFVB2","MCMC", "INLA","STAN"), 
       col=c("orange","blue","red", "purple","black"))

################################################################################
# Plot one of the Eta_10_Prediction
################################################################################
#MCMC
plot(density(eta_testdata_MCMC[10,]), main = "Density for Linear Predictor Prediction of test_Eta_10", col="red", ylim=c(0,2))
abline(v=eta_testdata_real[10,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_pred_mean[10] , sd = field_pred_sd[10])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

# Sampling INFVBphi 
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[10], sd=sqrt(cov_sample_eta_testdata_VB_all[10,10]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="orange", lty=1)

# Sampling INFVBssqphi 
ySeq_VBfixssqphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VBssqphi[10], sd=sqrt(cov_sample_eta_testdata_VBssqphi_all[10,10]))
lines(x=xSeq , y = ySeq_VBfixssqphi_sampling , typ="l" , col="blue", lty=1)


# Stan
lines(density(eta_testdata_MCMC_Stan[10,]))

legend("topright", lty=c(1,1,1), legend = c("INFVB1","INFVB2","MCMC", "INLA","STAN"), 
       col=c("orange","blue","red", "purple","black"))

################################################################################
# Plot one of the Eta_11_Prediction
################################################################################
#MCMC
plot(density(eta_testdata_MCMC[11,]), main = "Density for Linear Predictor Prediction of test_Eta_11", col="red", ylim=c(0,1.8))
abline(v=eta_testdata_real[11,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_pred_mean[11] , sd = field_pred_sd[11])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

# Sampling INFVBphi 
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[11], sd=sqrt(cov_sample_eta_testdata_VB_all[11,11]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="orange", lty=1)

# Sampling INFVBssqphi 
ySeq_VBfixssqphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VBssqphi[11], sd=sqrt(cov_sample_eta_testdata_VBssqphi_all[11,11]))
lines(x=xSeq , y = ySeq_VBfixssqphi_sampling , typ="l" , col="blue", lty=1)

legend("topright", lty=c(1,1,1), legend = c("INFVB1","INFVB2","MCMC", "INLA"), 
       col=c("orange","blue","red", "purple"))


# Stan
lines(density(eta_testdata_MCMC_Stan[11,]))

legend("topright", lty=c(1,1,1), legend = c("INFVB1","INFVB2","MCMC", "INLA","STAN"), 
       col=c("orange","blue","red", "purple","black"))


################################################################################
# Plot one of the Eta_20_Prediction
################################################################################
#MCMC
plot(density(eta_testdata_MCMC[20,]), main = "Density for Linear Predictor Prediction of test_Eta_20", col="red", ylim=c(0,2))
abline(v=eta_testdata_real[20,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_pred_mean[20] , sd = field_pred_sd[20])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

# Sampling INFVBphi 
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[20], sd=sqrt(cov_sample_eta_testdata_VB_all[20,20]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="orange", lty=1)

# Sampling INFVBssqphi 
ySeq_VBfixssqphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VBssqphi[20], sd=sqrt(cov_sample_eta_testdata_VBssqphi_all[20,20]))
lines(x=xSeq , y = ySeq_VBfixssqphi_sampling , typ="l" , col="blue", lty=1)


# Stan
lines(density(eta_testdata_MCMC_Stan[20,]))

legend("topright", lty=c(1,1,1), legend = c("INFVB1","INFVB2","MCMC", "INLA","STAN"), 
       col=c("orange","blue","red", "purple","black"))



aucVal
aucVal_VB
aucVal_INLA
aucVal_Stan

ptFinal
VB_ptFinal
comptTime_INLA[4]
totTime


aucVal
aucVal_INLA