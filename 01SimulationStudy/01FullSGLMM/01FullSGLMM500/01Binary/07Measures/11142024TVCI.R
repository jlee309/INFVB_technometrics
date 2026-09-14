################################################################################
################################################################################
# True Value in Credible Interval 
################################################################################
################################################################################
# 11142024 
################################################################################
################################################################################
# The real Data
getwd()
rm(list=ls())
load(file="../01DataGeneration/B_phi1_Gen500_1.RData")

eta_trainingdata_real <-X_trainingdata%*%beta+w_trainingdata
eta_testdata_real <-X_testdata%*%beta+w_testdata
################################################################################
################################################################################
# (1) MCMC methods 
load(file="../02MCMC/B_phi1_Gen500_1_MCMC.RData")
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
load(file="../06INLA/B_phi1_Gen500_1_INLA.RData")
library(INLA)
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
  lower_bound <- qnorm(0.025,mean = field_est_mean[k] , sd = field_est_sd[k])
  upper_bound <- qnorm(0.975,mean = field_est_mean[k] , sd = field_est_sd[k])
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
load(file="../03INFVBfixphi_ParallelComputing/B_phi1_Gen500_1_VBphiParallel.RData")
load(file="../08INFVBfixphi_SA/SA_B_phi1_Gen500_1_Full_INFVBphi.RData")

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
load(file="../09INFVBfixssqphi_SA/SA_B_phi1_Gen500_1_Full_INFVBssqphi.RData")

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

