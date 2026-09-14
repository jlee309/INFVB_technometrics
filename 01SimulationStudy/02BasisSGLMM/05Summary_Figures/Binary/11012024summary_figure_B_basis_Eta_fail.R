

#should be fixed  cov_sample_eta_testdata_VB_all!!!!!

setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/02BasisSGLMM/05Summary_Figures")
getwd()
rm(list=ls())
dev.off()
library(mvtnorm)
library(cli);library(viridis);library(fields);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm)
library(invgamma);library(pROC);library(Matrix);library(emulator);library(ggplot2)


# load(file="../../01DataGeneration/01phi01SpatialData20k.RData")
# load(file="../../01DataGeneration/01phi01EigenBasis20k.RData")
# load(file="../../02MCMC/Binary/01phi01Results_MCMC_Binary_Eigen.RData")
# load(file="../../03INFVB_ParallelComputing/Binary/01phi01ParallelINFVB_Binary_Eigen.RData")
# load(file="../../03INFVB_ParallelComputing/Binary/B_Basis_SA_INFVBssq.RData")
# load(file="../../04HybridMFVB/Binary/Basis_MFVB_SA.RData")
# load(file="../../06INLA/Binary/B_phi1_Gen500_1_INLA.RData")


load(file="../01DataGeneration/phi1_Gen25k_1SpatialData.RData")
load(file="../01DataGeneration/phi1_Gen25k_1EigenBasis.RData")
load(file="../02MCMC/Binary/Basis_B_phi1_Gen25k_1_MCMC.RData")
load(file="../03INFVB_ParallelComputing/Binary/Basis_B_phi1_Gen25k_1_VBssqParallel.RData")
load(file="../07INFVB_SA/Binary/SA_Basis_B_phi1_Gen25k_1_VBssqParallel.RData")
# load(file="../04HybridMFVB/Binary/Basis_B_phi1_Gen25k_1_HMFVB.RData")






iter.thetaphi<-1000


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
plot(density(eta_trainingdata_MCMC[3,]), main = "Density for Linear Predictor Eta_3", col="red")
abline(v=eta_trainingdata_real[3,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
# ySeq<-dnorm(x = xSeq , mean = field_est_mean[3] , sd = field_est_sd[3])
# lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VB[3], sd=sqrt(cov_sample_eta_trainingdata_VB_all[3,3]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Est[3,], sd=sqrt(cov_sample_eta_trainingdata_MFVB[3,3]))
lines(x=xSeq , y = ySeq_MFVB_sampling , typ="l" , col="orange", lty=1)

legend("topright", lty=c(1,1,1), legend = c("SA_VBfix1","MCMC", "INLA","SA_MFVB"), 
       col=c("blue","red", "purple","orange"))


################################################################################
# Plot one of the Eta_10_estimation
################################################################################
# Eta_10
#MCMC
plot(density(eta_trainingdata_MCMC[10,]), main = "Density for Linear Predictor Eta_10", col="red")
abline(v=eta_trainingdata_real[10,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_est_mean[10] , sd = field_est_sd[10])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VB[10], sd=sqrt(cov_sample_eta_trainingdata_VB_all[10,10]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Est[10,], sd=sqrt(cov_sample_eta_trainingdata_MFVB[10,10]))
lines(x=xSeq , y = ySeq_MFVB_sampling , typ="l" , col="orange", lty=1)

legend("topright", lty=c(1,1,1), legend = c("SA_VBfix1","MCMC", "INLA","SA_MFVB"), 
       col=c("blue","red", "purple","orange"))
################################################################################
# Plot one of the Eta_25_estimation
################################################################################
# Eta_25
#MCMC
plot(density(eta_trainingdata_MCMC[25,]), main = "Density for Linear Predictor Eta_25", col="red")
abline(v=eta_trainingdata_real[25,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_est_mean[25] , sd = field_est_sd[25])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VB[25], sd=sqrt(cov_sample_eta_trainingdata_VB_all[25,25]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Est[25,], sd=sqrt(cov_sample_eta_trainingdata_MFVB[25,25]))
lines(x=xSeq , y = ySeq_MFVB_sampling , typ="l" , col="orange", lty=1)

legend("topright", lty=c(1,1,1), legend = c("SA_VBfix1","MCMC", "INLA","SA_MFVB"), 
       col=c("blue","red", "purple","orange"))

################################################################################
# Plot one of the Eta_40_estimation
################################################################################
# Eta_40
#MCMC
plot(density(eta_trainingdata_MCMC[40,]), main = "Density for Linear Predictor Eta_40", col="red")
abline(v=eta_trainingdata_real[40,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_est_mean[40] , sd = field_est_sd[40])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VB[40], sd=sqrt(cov_sample_eta_trainingdata_VB_all[40,40]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Est[40,], sd=sqrt(cov_sample_eta_trainingdata_MFVB[40,40]))
lines(x=xSeq , y = ySeq_MFVB_sampling , typ="l" , col="orange", lty=1)

legend("topright", lty=c(1,1,1), legend = c("SA_VBfix1","MCMC", "INLA","SA_MFVB"), 
       col=c("blue","red", "purple","orange"))















################################################################################
# (5-1) Plotting Etas in the test datasets
################################################################################
################################################################################
# Plot one of the Eta_11_Prediction
################################################################################
#MCMC
plot(density(eta_testdata_MCMC[3,]), main = "Density for Linear Predictor Prediction of test_Eta_3", col="red")
abline(v=eta_testdata_real[3,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_pred_mean[3] , sd = field_pred_sd[3])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[3], sd=sqrt(cov_sample_eta_testdata_VB_all[3,3]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Pred[3,], sd=sqrt(cov_sample_eta_testdata_MFVB[3,3]))
lines(x=xSeq , y = ySeq_MFVB_sampling , typ="l" , col="orange", lty=1)

legend("topright", lty=c(1,1,1), legend = c("SA_VBfix1","MCMC", "INLA","SA_MFVB"), 
       col=c("blue","red", "purple","orange"))
################################################################################
# Plot one of the Eta_10_Prediction
################################################################################
#MCMC
plot(density(eta_testdata_MCMC[10,]), main = "Density for Linear Predictor Prediction of test_Eta_10", col="red")
abline(v=eta_testdata_real[10,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_pred_mean[10] , sd = field_pred_sd[10])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[10], sd=sqrt(cov_sample_eta_testdata_VB_all[10,10]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Pred[10,], sd=sqrt(cov_sample_eta_testdata_MFVB[10,10]))
lines(x=xSeq , y = ySeq_MFVB_sampling , typ="l" , col="orange", lty=1)

legend("topright", lty=c(1,1,1), legend = c("SA_VBfix1","MCMC", "INLA","SA_MFVB"), 
       col=c("blue","red", "purple","orange"))

################################################################################
# Plot one of the Eta_11_Prediction
################################################################################
#MCMC
plot(density(eta_testdata_MCMC[11,]), main = "Density for Linear Predictor Prediction of test_Eta_11", col="red")
abline(v=eta_testdata_real[11,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_pred_mean[11] , sd = field_pred_sd[11])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[11], sd=sqrt(cov_sample_eta_testdata_VB_all[11,11]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Pred[11,], sd=sqrt(cov_sample_eta_testdata_MFVB[11,11]))
lines(x=xSeq , y = ySeq_MFVB_sampling , typ="l" , col="orange", lty=1)

legend("topright", lty=c(1,1,1), legend = c("SA_VBfix1","MCMC", "INLA","SA_MFVB"), 
       col=c("blue","red", "purple","orange"))

################################################################################
# Plot one of the Eta_20_Prediction
################################################################################
#MCMC
plot(density(eta_testdata_MCMC[20,]), main = "Density for Linear Predictor Prediction of test_Eta_20", col="red")
abline(v=eta_testdata_real[20,], col="green")

# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = field_pred_mean[20] , sd = field_pred_sd[20])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)


#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[20], sd=sqrt(cov_sample_eta_testdata_VB_all[20,20]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Pred[20,], sd=sqrt(cov_sample_eta_testdata_MFVB[20,20]))
lines(x=xSeq , y = ySeq_MFVB_sampling , typ="l" , col="orange", lty=1)

legend("topright", lty=c(1,1,1), legend = c("SA_VBfix1","MCMC", "INLA","SA_MFVB"), 
       col=c("blue","red", "purple","orange"))

