# problems TVCI in VB and MCMC
# Eta_3
#MCMC
plot(density(eta_trainingdata_MCMC[3,]), main = "Density for Linear Predictor Eta_3", col="red", xlim=c(-2,2))
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



# Eta_1
#MCMC
plot(density(eta_trainingdata_MCMC[1,]), main = "Density for Linear Predictor Eta_1", col="red", xlim=c(-2,2))
abline(v=eta_trainingdata_real[1,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = linear_predictor_mean[1] , sd = linear_predictor_sd[1])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_est_VB[1], sd=sqrt(cov_sample_eta_trainingdata_VB_all[1,1]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Est[1], sd=sqrt(cov_sample_eta_trainingdata_MFVB[1,1]))
lines(x=xSeq , y = ySeq_MFVB_sampling , typ="l" , col="orange", lty=1)

legend("topright", lty=c(1,1,1), legend = c("SA_VBfix1","MCMC", "INLA","SA_MFVB"), 
       col=c("blue","red", "purple","orange"))

obs_RMSPE
obs_RMSPE_INLA



################################################################################
# Plot one of the Eta_11_Prediction
################################################################################
#MCMC
plot(density(eta_testdata_MCMC[11,]), main = "Density for Linear Predictor Prediction of test_Eta_11", col="red", xlim=c(-2,2))
abline(v=eta_testdata_real[11,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = linear_predictor_mean[(0.8*n)+11] , sd = linear_predictor_sd[(0.8*n)+11])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[11], sd=sqrt(cov_sample_eta_testdata_VB_all[11,11]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Pred[11], sd=sqrt(cov_sample_eta_testdata_MFVB[11,11]))
lines(x=xSeq , y = ySeq_MFVB_sampling , typ="l" , col="orange", lty=1)

legend("topright", lty=c(1,1,1), legend = c("SA_VBfix1","MCMC", "INLA","SA_MFVB"), 
       col=c("blue","red", "purple","orange"))

obs_RMSPE
obs_RMSPE_INLA
obs_RMSPE_VB

TotTime_INLA
