
rm(list=ls())
setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/02BasisSGLMM/05Summary_Figures/Poisson")
dev.off()
getwd()

library(cli);library(viridis);library(fields);library(pgdraw);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm);
library(pROC);library(Matrix);library(emulator);

################################################################################
################################################################################
################################################################################
################################################################################

# Basis_C_phi1_Gen25k_1_MCMC.RData
# Basis_C_phi1_Gen25k_1_VBssqParallel.RData
# SA_Basis_C_phi2_Gen25k_50_HMFVB.RData
# Basis_C_phi1_Gen25k_4_HMFVB.RData

load(file="../../01DataGeneration/phi1_Gen25k_1SpatialData.RData")
load(file="../../01DataGeneration/phi1_Gen25k_1EigenBasis.RData")
load(file="../../02MCMC/Poisson/Basis_C_phi1_Gen25k_1_MCMC.RData")
load(file="../../03INFVB_ParallelComputing/Poisson/Basis_C_phi1_Gen25k_1_VBssqParallel.RData")
load(file="../../04HybridMFVB/Poisson/Basis_C_phi1_Gen25k_1_HMFVB.RData")
load(file="../../07INFVB_SA/Poisson/SA_Basis_C_phi1_Gen25k_1_VBssqParallel.RData") 
load(file="../../08HybridMFVB_SA/Poisson/SA_Basis_C_phi1_Gen25k_1_HMFVB.RData") 
load(file="../../06INLA/Poisson/Basis_C_phi1_Gen25k_1_INLA.RData")

################################################################################
# (5) Plotting Etas in the training datasets
################################################################################
# The real Data
eta_trainingdata_real <-X_trainingdata%*%beta+w_trainingdata
eta_testdata_real <-X_testdata%*%beta+w_testdata









################################################################################
################################################################################
## Comparison using Plot
# Setting for overall MFVB
j<-nrow(MFVB_Meangamma)
p<-2
burninperiod<-ceiling(n.iteration*0.1)
################################################################################
################################################################################
# ggplot
################################################################################
################################################################################
################################################################################
################################################################################
# Beta 1
################################################################################
library(ggplot2)
MFVB_beta1<-dnorm(x=xSeqb1, mean = MFVB_Meangamma[j,1] ,sd =  sqrt(MFVB_covMatgamma[j,1]))
df_beta1_1<-data.frame('beta1mcmc'=parameterMatrix[-(1:burninperiod),1])
df_beta1_2<-data.frame(newdensityBeta1,xSeqb1)
df_beta1_3<-data.frame(MFVB_beta1,xSeqb1)

ggplot(df_beta1_1,aes(x=beta1mcmc))+
  geom_histogram(aes(y=..density..,color="MCMC"),binwidth =0.005)+
  geom_line(data=df_beta1_2, aes(x=xSeqb1,y=newdensityBeta1,color="INFVB"),size=1.5)+
  geom_line(data=df_beta1_3, aes(x=xSeqb1,y=MFVB_beta1,color="MFVB"),size=1.5,linetype=5)+
  scale_color_manual(name = "",values = c(MCMC = "black", INFVB = "blue",MFVB="red"))+
  coord_cartesian(xlim=c(0.8,1.2))+
  geom_vline(xintercept =beta[1], color="green",size=1.5)+
  theme(legend.position = "right")+
  ggtitle("Beta1")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=12,
                                 color="black",
                                 angle=0),
        axis.text.x=element_text(size=12,
                                 color="black",
                                 angle=0),
        plot.title = element_text(size=20,
                                  hjust = 0.5))

################################################################################
## Beta 2
################################################################################
MFVB_beta2<-dnorm(x=xSeqb2, mean = MFVB_Meangamma[j,2] ,sd =  sqrt(MFVB_covMatgamma[j,2]))
df_beta2_1<-data.frame('beta2mcmc'=parameterMatrix[-(1:burninperiod),2])
df_beta2_2<-data.frame(newdensityBeta2,xSeqb2)
df_beta2_3<-data.frame(MFVB_beta2,xSeqb2)

ggplot(df_beta2_1,aes(x=beta2mcmc))+
  geom_histogram(aes(y=..density..,color="MCMC"),binwidth =0.005)+
  geom_line(data=df_beta2_2, aes(x=xSeqb1,y=newdensityBeta2,color="INFVB"),size=1.5)+
  geom_line(data=df_beta2_3, aes(x=xSeqb1,y=MFVB_beta2,color="MFVB"),size=1.5,linetype=5)+
  scale_color_manual(name = "",values = c(MCMC = "black", INFVB = "blue",MFVB="red"))+
  coord_cartesian(xlim=c(0.8,1.2))+
  geom_vline(xintercept =beta[2], color="green",size=1.5)+
  theme(legend.position = "right")+
  ggtitle("Beta2")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=12,
                                 color="black",
                                 angle=0),
        axis.text.x=element_text(size=12,
                                 color="black",
                                 angle=0),
        plot.title = element_text(size=20,
                                  hjust = 0.5))


################################################################################
## Sigma2
################################################################################
xSeqsigma2<-seq(0,max(x_sigma2), length.out=iter.thetaphi)
set.seed(5153527)
foosigma2<-sample(x_sigma2,size = iter.thetaphi^2, replace = T, prob = ELBOvectorweights)

MFVB_sigma2<-dinvgamma(x=xSeqsigma2, shape = MFVB_iter[j,1] ,rate=MFVB_iter[j,2])
df_sigma2_1<-data.frame('sigma2mcmc'=parameterMatrix[-(1:burninperiod),3])
df_sigma2_2<-data.frame(foosigma2)
df_sigma2_3<-data.frame(MFVB_sigma2,xSeqsigma2)

ggplot(df_sigma2_1,aes(x=sigma2mcmc))+
  geom_histogram(aes(y=..density..,color="MCMC"),binwidth =20)+
  geom_density(data=df_sigma2_2, aes(x=foosigma2, color="INFVB"),size=1.5)+
  geom_line(data=df_sigma2_3, aes(x=xSeqsigma2,y=MFVB_sigma2,color="MFVB"),size=1.5,linetype=5)+
  scale_color_manual(name = "",values = c(MCMC = "black", INFVB = "blue",MFVB="red"))+
  #geom_vline(xintercept =sigma2, color="green",size=1)+
  coord_cartesian(xlim=c(0,1000))+
  theme(legend.position = "right")+
  ggtitle("Sigma2")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=12,
                                 color="black",
                                 angle=0),
        axis.text.x=element_text(size=12,
                                 color="black",
                                 angle=0),
        plot.title = element_text(size=20,
                                  hjust = 0.5))

################################################################################
## Delta3
################################################################################
MFVB_delta3<-dnorm(x=xSeqdelta3, mean=MFVB_Meangamma[j,p+3], sd=sqrt(MFVB_covMatgamma[j,p+3]))
df_delta3_1<-data.frame('delta3mcmc'=deltaMatrix[-(1:burninperiod),3])
df_delta3_2<-data.frame(newdensitydelta3,xSeqdelta3)
df_delta3_3<-data.frame(MFVB_delta3,xSeqdelta3)

ggplot(df_delta3_1,aes(x=delta3mcmc))+
  geom_histogram(aes(y=..density..,color="MCMC"),binwidth =0.2)+
  geom_line(data=df_delta3_2, aes(x=xSeqdelta3,y=newdensitydelta3,color="INFVB"),size=1.5)+
  geom_line(data=df_delta3_3, aes(x=xSeqdelta3,y=MFVB_delta3,color="MFVB"),size=1.5,linetype=5)+
  scale_color_manual(name = "",values = c(MCMC = "black", INFVB = "blue",MFVB="red"))+
  coord_cartesian(xlim=c(-10,5))+
  #geom_vline(xintercept =beta[2], color="green",size=1.5)+
  theme(legend.position = "right")+
  ggtitle("delta3")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=12,
                                 color="black",
                                 angle=0),
        axis.text.x=element_text(size=12,
                                 color="black",
                                 angle=0),
        plot.title = element_text(size=20,
                                  hjust = 0.5))
################################################################################
################################################################################
################################################################################
################################################################################
#loglambda_testdata-Truth
df_loglambda_testdata<-data.frame(log(lambda_testdata),locations_testdata)
head(df_loglambda_testdata)

ggplot(df_loglambda_testdata)+
  geom_point(aes(x=df_loglambda_testdata[,2] , y=df_loglambda_testdata[,3] , color=log(lambda_testdata)), shape=7)+
  scale_color_gradientn(colours = rainbow(3))+
  # labs(title="Probability - Truth", y="", x="")+
  theme(legend.position = "right",
        legend.title=element_blank() ) +
  theme(axis.title.x=element_text(face="plain",
                                  size=15,
                                  hjust=0.5,
                                  color="black"),
        axis.title.y=element_text(face="plain",
                                  size=15,
                                  hjust=0.5,
                                  color="black"),
        axis.text.y=element_text(size=12,
                                 color="black",
                                 angle=0),
        axis.text.x=element_text(size=12,
                                 color="black",
                                 angle=0),
        plot.title = element_text(size=20,
                                  hjust = 0.5))

################################################################################
################################################################################
#log(lambda_prediction)-MCMC
df_loglambda_prediction<-data.frame(log(lambda_prediction),locations_testdata)
head(df_loglambda_prediction)

ggplot(df_loglambda_prediction)+
  geom_point(aes(x=df_loglambda_prediction[,2] , y=df_loglambda_prediction[,3] , color=log(lambda_prediction)), shape=7)+
  scale_color_gradientn(colours = rainbow(3))+
  # labs(title="Probability - MCMC", y="", x="")+
  theme(legend.position = "right",
        legend.title=element_blank() ) +
  theme(axis.title.x=element_text(face="plain",
                                  size=15,
                                  hjust=0.5,
                                  color="black"),
        axis.title.y=element_text(face="plain",
                                  size=15,
                                  hjust=0.5,
                                  color="black"),
        axis.text.y=element_text(size=12,
                                 color="black",
                                 angle=0),
        axis.text.x=element_text(size=12,
                                 color="black",
                                 angle=0),
        plot.title = element_text(size=20,
                                  hjust = 0.5))

################################################################################
################################################################################
#prob_testdata-INFVB
df_loglambda_prediction_VB<-data.frame(log(lambda_prediction_VB),locations_testdata)
head(df_loglambda_prediction_VB)

ggplot(df_loglambda_prediction_VB)+
  geom_point(aes(x=df_loglambda_prediction_VB[,2] , y=df_loglambda_prediction_VB[,3] , color=log(lambda_prediction_VB)), shape=7)+
  scale_color_gradientn(colours = rainbow(3))+
  #labs(title="Probability - VB", y="", x="")+
  theme(legend.position = "right",
        legend.title=element_blank() ) +
  theme(axis.title.x=element_text(face="plain",
                                  size=15,
                                  hjust=0.5,
                                  color="black"),
        axis.title.y=element_text(face="plain",
                                  size=15,
                                  hjust=0.5,
                                  color="black"),
        axis.text.y=element_text(size=12,
                                 color="black",
                                 angle=0),
        axis.text.x=element_text(size=12,
                                 color="black",
                                 angle=0),
        plot.title = element_text(size=20,
                                  hjust = 0.5))
################################################################################
################################################################################
df_loglambda_prediction_MFVB<-data.frame(log(lambda_prediction_MFVB),locations_testdata)
head(df_loglambda_prediction_MFVB)

ggplot(df_loglambda_prediction_MFVB)+
  geom_point(aes(x=df_loglambda_prediction_MFVB[,2] , y=df_loglambda_prediction_MFVB[,3] , color=log(lambda_prediction_MFVB)), shape=7)+
  scale_color_gradientn(colours = rainbow(3))+
  #labs(title="Probability - VB", y="", x="")+
  theme(legend.position = "right",
        legend.title=element_blank() ) +
  theme(axis.title.x=element_text(face="plain",
                                  size=15,
                                  hjust=0.5,
                                  color="black"),
        axis.title.y=element_text(face="plain",
                                  size=15,
                                  hjust=0.5,
                                  color="black"),
        axis.text.y=element_text(size=12,
                                 color="black",
                                 angle=0),
        axis.text.x=element_text(size=12,
                                 color="black",
                                 angle=0),
        plot.title = element_text(size=20,
                                  hjust = 0.5))




################################################################################
################################################################################
## RMSPE FINAL RESULT Comparison between MCMC vs VB 
################################################################################
################################################################################
obs_RMSPE # RMSPE - MCMC
obs_RMSPE_MFVB # RMSPE - MFVB
obs_RMSPE_VB # RMSPE - INFVB

ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)


################################################################################
################################################################################

################################################################################
# Plot one of the Eta_11_Prediction
################################################################################
# The real Data
eta_trainingdata_real <-X_trainingdata%*%beta+w_trainingdata
eta_testdata_real <-X_testdata%*%beta+w_testdata

#MCMC
plot(density(eta_testdata_MCMC[3,]), main = "Density for Linear Predictor Prediction of test_Eta_3", col="red", xlim=c(-1,2))
abline(v=eta_testdata_real[3,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = linear_predictor_mean[3] , sd = linear_predictor_sd[3])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[3], sd=sqrt(cov_sample_eta_testdata_VB_all[3,3]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Pred[3], sd=sqrt(cov_sample_eta_testdata_MFVB[3,3]))
lines(x=xSeq , y = ySeq_MFVB_sampling , typ="l" , col="orange", lty=1)

legend("topright", lty=c(1,1,1), legend = c("SA_VBfix1","MCMC", "INLA","SA_MFVB"), 
       col=c("blue","red", "purple","orange"))




#MCMC
plot(density(eta_testdata_MCMC[11,]), main = "Density for Linear Predictor Prediction of test_Eta_11", col="red", xlim=c(-1,2))
abline(v=eta_testdata_real[11,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = linear_predictor_mean[11] , sd = linear_predictor_sd[11])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[11], sd=sqrt(cov_sample_eta_testdata_VB_all[11,11]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-11,11,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Pred[11], sd=sqrt(cov_sample_eta_testdata_MFVB[11,11]))
lines(x=xSeq , y = ySeq_MFVB_sampling , typ="l" , col="orange", lty=1)

legend("topright", lty=c(1,1,1), legend = c("SA_VBfix1","MCMC", "INLA","SA_MFVB"), 
       col=c("blue","red", "purple","orange"))



#MCMC
plot(density(exp(eta_testdata_MCMC[11,])), main = "Density for Linear Predictor Prediction of test_Eta_11", col="red")
abline(v=eta_testdata_real[11,], col="green")
# INLA
xSeq<-seq(-3,3,length.out=iter.thetaphi)
ySeq<-dnorm(x = xSeq , mean = linear_predictor_mean[11] , sd = linear_predictor_sd[11])
lines(x=xSeq , y = ySeq , typ="l" , col="purple", lty=1)

#INFVB(ssq)_Sampling
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[11], sd=sqrt(cov_sample_eta_testdata_VB_all[11,11]))
lines(x=xSeq , y = ySeq_VBfixphi_sampling , typ="l" , col="blue", lty=1)

# HMFVB
xSeq<-seq(-11,11,length.out=iter.thetaphi)
ySeq_MFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Pred[11], sd=sqrt(cov_sample_eta_testdata_MFVB[11,11]))
lines(x=xSeq , y = ySeq_MFVB_sampling , typ="l" , col="orange", lty=1)

legend("topright", lty=c(1,1,1), legend = c("SA_VBfix1","MCMC", "INLA","SA_MFVB"), 
       col=c("blue","red", "purple","orange"))


