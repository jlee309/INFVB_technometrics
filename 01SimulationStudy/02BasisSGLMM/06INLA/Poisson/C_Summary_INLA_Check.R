


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
dev.off()
rm(list=ls())
file_name1 <- paste0("../../01DataGeneration/phi1_Gen25k_1EigenBasis.RData")
file_name2 <-paste0("../../01DataGeneration/phi1_Gen25k_1SpatialData.RData")
load(file_name1)
load(file_name2)
load(file="../../02MCMC/Poisson/Basis_C_phi1_Gen25k_1_MCMC.RData")
load(file="../../06INLA/Poisson/Basis_C_phi1_Gen25k_1_INLA.RData")
load(file="../../03INFVB_ParallelComputing/Poisson/Basis_C_phi1_Gen25k_1_VBssqParallel.RData")
load(file="../../07INFVB_SA/Poisson/SA_Basis_C_phi1_Gen25k_1_VBssqParallel.RData")
# load(file="../../08HybridMFVB_SA/Poisson/SA_Basis_C_phi1_Gen25k_1_HMFVB.RData")
# load(file="../../06INLA/Poisson/B_phi1_Gen500_1_INLA.RData")
# load(file="../../10Stan/Poisson/mcmc_Stan_small.RData")

################################################################################################
#Check INLA with MCMC
# Create the first density plot
burnin<-10000
plot(density(parameterMatrix[-burnin, 1]), col = "blue", main = "Density Plots", xlab = "Parameter Value", ylab = "Density")
# INLA
xSeq<-seq(-3,3,length.out=1000)
ySeq<-dnorm(x = xSeq , mean = INLA_Basis_Mean[1] , sd = INLA_Basis_SD[1]) 
lines(x=xSeq , y = ySeq , typ="l" , col="red", lty=1)
# Optionally, add a legend
legend("topright", legend = c("MCMC", "INLA"), col = c("blue", "red"), lty = 1)
################################################################################################
#Check INLA with MCMC
# Create the first density plot
plot(density(parameterMatrix[-burnin, 2]), col = "blue", main = "Density Plots", xlab = "Parameter Value", ylab = "Density")
# INLA
xSeq<-seq(-3,3,length.out=1000)
ySeq<-dnorm(x = xSeq , mean = INLA_Basis_Mean[2] , sd = INLA_Basis_SD[2]) 
lines(x=xSeq , y = ySeq , typ="l" , col="red", lty=1)
# Optionally, add a legend
legend("topright", legend = c("MCMC", "INLA"), col = c("blue", "red"), lty = 1)
################################################################################################
#Delta 
plot(density(deltaMatrix[-burnin,2]), col = "blue", main = "Density Plots", xlab = "Parameter Value", ylab = "Density")
# INLA
xSeq<-seq(-10,10,length.out=1000)
ySeq<-dnorm(x = xSeq , mean = INLA_Basis_Mean[(4)] , sd = INLA_Basis_SD[(4)]) 
lines(x=xSeq , y = ySeq , typ="l" , col="red", lty=1)
# Optionally, add a legend
legend("topright", legend = c("MCMC", "INLA"), col = c("blue", "red"), lty = 1)



#Delta 
plot(density(deltaMatrix[-burnin,8]), col = "blue", main = "Density Plots", xlab = "Parameter Value", ylab = "Density")
# INLA
xSeq<-seq(-10,10,length.out=1000)
ySeq<-dnorm(x = xSeq , mean = INLA_Basis_Mean[(p+8)] , sd = INLA_Basis_SD[(p+8)]) 
lines(x=xSeq , y = ySeq , typ="l" , col="red", lty=1)
# Optionally, add a legend
legend("topright", legend = c("MCMC", "INLA"), col = c("blue", "red"), lty = 1)




#Delta 
plot(density(deltaMatrix[-burnin,10]), col = "blue", main = "Density Plots", xlab = "Parameter Value", ylab = "Density")
# INLA
xSeq<-seq(-10,10,length.out=1000)
ySeq<-dnorm(x = xSeq , mean = INLA_Basis_Mean[(p+10)] , sd = INLA_Basis_SD[(p+10)]) 
lines(x=xSeq , y = ySeq , typ="l" , col="red", lty=1)
# Optionally, add a legend
legend("topright", legend = c("MCMC", "INLA"), col = c("blue", "red"), lty = 1)





################################################################################
################################################################################
################################################################################
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
ySeq<-dnorm(x = xSeq , mean = linear_predictor_mean[3] , sd = linear_predictor_sd[3])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VB[3], sd=sqrt(cov_sample_eta_trainingdata_VB_all[3,3]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Est[3], sd=sqrt(cov_sample_eta_trainingdata_MFVB[3,3]))
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
ySeq<-dnorm(x = xSeq , mean = linear_predictor_mean[10] , sd = linear_predictor_sd[10])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VB[10], sd=sqrt(cov_sample_eta_trainingdata_VB_all[10,10]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Est[10], sd=sqrt(cov_sample_eta_trainingdata_MFVB[10,10]))
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
plot(density(eta_testdata_MCMC[3,]), main = "Density for Linear Predictor Prediction of test_Eta_3", col="red", ylim=c(0,2))
abline(v=eta_testdata_real[3,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = linear_predictor_mean[403] , sd = linear_predictor_sd[403])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)
# plot(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[3], sd=sqrt(cov_sample_eta_testdata_VB_all[3,3]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Pred[3], sd=sqrt(cov_sample_eta_testdata_MFVB[3,3]))
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
rth<-50
plot(density(eta_testdata_MCMC[rth,]), main = "Density for Linear Predictor Prediction of test_Eta_3", col="red", xlim=c(0-10,10))
abline(v=eta_testdata_real[rth,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = linear_predictor_mean[(0.8*n)+rth] , sd = linear_predictor_sd[(0.8*n)+rth])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)
# plot(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[rth], sd=sqrt(cov_sample_eta_testdata_VB_all[rth,rth]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

