################################################################################
################################################################################
# 11062025 Summary for Binary 50 Basis functions
################################################################################
################################################################################
rm(list=ls())
setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/02BasisSGLMM/05Summary_Figures/Binary")
# dev.off()
getwd()

library(cli);library(viridis);library(fields);library(pgdraw);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm);
library(pROC);library(Matrix);library(emulator);

################################################################################
################################################################################
################################################################################
################################################################################
# dev.off()
# load(file="../../01DataGeneration/phi1_Gen25k_1SpatialData.RData")
# load(file="../../01DataGeneration/phi1_Gen25k_1EigenBasis.RData")
# load(file="../../02MCMC/Binary/50Basis_B_phi1_Gen25k_1_MCMC.RData")
# load(file="../../03INFVB_ParallelComputing/Binary/50Basis_B_phi1_Gen25k_1_VBssqParallel.RData")
# load(file="../../04HybridMFVB/Binary/50Basis_B_phi1_Gen25k_1_HMFVB.RData")
# load(file="../../06INLA/Binary/50Basis_B_phi1_Gen25k_1_INLA.RData")
# load(file="../../07INFVB_SA/Binary/SA_50Basis_B_phi1_Gen25k_1_VBssqParallel.RData") 
# load(file="../../08HybridMFVB_SA/Binary/SA_50Basis_B_phi1_Gen25k_1_HMFVB.RData") 
# load(file="../../10Stan/Binary/50Basis_B_phi1_Gen25k_1_MCMC_Stan.RData") 


# load(file="../../01DataGeneration/phi2_Gen25k_1SpatialData.RData")
# load(file="../../01DataGeneration/phi2_Gen25k_1EigenBasis.RData")
# load(file="../../02MCMC/Binary/50Basis_B_phi2_Gen25k_1_MCMC.RData")
# load(file="../../03INFVB_ParallelComputing/Binary/50Basis_B_phi2_Gen25k_1_VBssqParallel.RData")
# load(file="../../04HybridMFVB/Binary/50Basis_B_phi2_Gen25k_1_HMFVB.RData")
# load(file="../../06INLA/Binary/50Basis_B_phi2_Gen25k_1_INLA.RData")
# load(file="../../07INFVB_SA/Binary/SA_50Basis_B_phi2_Gen25k_1_VBssqParallel.RData") 
# load(file="../../08HybridMFVB_SA/Binary/SA_50Basis_B_phi2_Gen25k_1_HMFVB.RData") 
# load(file="../../10Stan/Binary/50Basis_B_phi2_Gen25k_1_MCMC_Stan.RData") 



load(file="../../01DataGeneration/phi3_Gen25k_1SpatialData.RData")
load(file="../../01DataGeneration/phi3_Gen25k_1EigenBasis.RData")
load(file="../../02MCMC/Binary/50Basis_B_phi3_Gen25k_1_MCMC.RData")
load(file="../../03INFVB_ParallelComputing/Binary/50Basis_B_phi3_Gen25k_1_VBssqParallel.RData")
load(file="../../04HybridMFVB/Binary/50Basis_B_phi3_Gen25k_1_HMFVB.RData")
load(file="../../06INLA/Binary/50Basis_B_phi3_Gen25k_1_INLA.RData")
load(file="../../07INFVB_SA/Binary/SA_50Basis_B_phi3_Gen25k_1_VBssqParallel.RData") 
load(file="../../08HybridMFVB_SA/Binary/SA_50Basis_B_phi3_Gen25k_1_HMFVB.RData") 
# load(file="../../10Stan/Binary/50Basis_B_phi3_Gen25k_1_MCMC_Stan.RData") 


# load(file="../../01DataGeneration/phi4_Gen25k_1SpatialData.RData")
# load(file="../../01DataGeneration/phi4_Gen25k_1EigenBasis.RData")
# load(file="../../02MCMC/Binary/50Basis_B_phi4_Gen25k_1_MCMC.RData")
# load(file="../../03INFVB_ParallelComputing/Binary/50Basis_B_phi4_Gen25k_1_VBssqParallel.RData")
# load(file="../../04HybridMFVB/Binary/50Basis_B_phi4_Gen25k_1_HMFVB.RData")
# load(file="../../06INLA/Binary/50Basis_B_phi4_Gen25k_1_INLA.RData")
# load(file="../../07INFVB_SA/Binary/SA_50Basis_B_phi4_Gen25k_1_VBssqParallel.RData") 
# load(file="../../08HybridMFVB_SA/Binary/SA_50Basis_B_phi4_Gen25k_1_HMFVB.RData") 
# load(file="../../10Stan/Binary/50Basis_B_phi4_Gen25k_1_MCMC_Stan.RData") 

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
j<-nrow(MFVB_MeanV)
# numofbasis<-20
p<-2
burninperiod<-ceiling(n.iteration*0.1)








# ggplot
################################################################################
################################################################################
################################################################################
################################################################################
# Beta 1
################################################################################
library(ggplot2)
MFVB_beta1<-dnorm(x=xSeqb1, mean = MFVB_MeanV[j,1] ,sd =  sqrt(MFVB_covMatV[j,1]))
df_beta1_1<-data.frame('beta1mcmc'=parameterMatrix[-(1:burninperiod),1])
df_beta1_2<-data.frame(newdensityBeta1,xSeqb1)
df_beta1_3<-data.frame(MFVB_beta1,xSeqb1)
# df_beta1_4<-data.frame('beta1stan'=parMat[,1])

# ggplot(df_beta1_1,aes(x=beta1mcmc))+
#   geom_histogram(aes(y=..density..,color="MCMC"),binwidth =0.005)+
#   geom_line(data=df_beta1_2, aes(x=xSeqb1,y=newdensityBeta1,color="INFVB"),size=1.5)+
#   geom_line(data=df_beta1_3, aes(x=xSeqb1,y=MFVB_beta1,color="MFVB"),size=1.5,linetype=5)+
#   scale_color_manual(name = "",values = c(MCMC = "black", INFVB = "blue",MFVB="red"))+
#   coord_cartesian(xlim=c(0.8,1.3))+
#   geom_vline(xintercept =beta[1], color="green",size=1.5)+
#   theme(legend.position = "right")+
#   ggtitle("Beta1")+
#   theme(axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_text(size=12,
#                                  color="black",
#                                  angle=0),
#         axis.text.x=element_text(size=12,
#                                  color="black",
#                                  angle=0),
#         plot.title = element_text(size=20,
#                                   hjust = 0.5))


ggplot(df_beta1_1,aes(x=beta1mcmc))+
  geom_histogram(aes(y=..density..,color="MCMC"),binwidth =0.005)+
  geom_line(data=df_beta1_2, aes(x=xSeqb1,y=newdensityBeta1,color="INFVB"),size=1.5)+
  geom_line(data=df_beta1_3, aes(x=xSeqb1,y=MFVB_beta1,color="MFVB"),size=1.5,linetype=5)+
  # geom_density(data = df_beta1_4, aes(x = beta1stan, fill = "STAN", color = "STAN"), alpha = 0.2, linetype = 1, size = 1.5) +
  scale_color_manual(name = "",values = c(MCMC = "black", INFVB = "blue",MFVB="red",STAN="yellow"))+
  coord_cartesian(xlim=c(0.85,1.1))+
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
MFVB_beta2<-dnorm(x=xSeqb2, mean = MFVB_MeanV[j,2] ,sd =  sqrt(MFVB_covMatV[j,2]))
df_beta2_1<-data.frame('beta2mcmc'=parameterMatrix[-(1:burninperiod),2])
df_beta2_2<-data.frame(newdensityBeta2,xSeqb2)
df_beta2_3<-data.frame(MFVB_beta2,xSeqb2)
# df_beta2_4<-data.frame('beta2stan'=parMat[,2])

# ggplot(df_beta2_1,aes(x=beta2mcmc))+
#   geom_histogram(aes(y=..density..,color="MCMC"),binwidth =0.005)+
#   geom_line(data=df_beta2_2, aes(x=xSeqb1,y=newdensityBeta2,color="INFVB"),size=1.5)+
#   geom_line(data=df_beta2_3, aes(x=xSeqb1,y=MFVB_beta2,color="MFVB"),size=1.5,linetype=5)+
#   scale_color_manual(name = "",values = c(MCMC = "black", INFVB = "blue",MFVB="red"))+
#   coord_cartesian(xlim=c(0.8,1.3))+
#   geom_vline(xintercept =beta[2], color="green",size=1.5)+
#   theme(legend.position = "right")+
#   ggtitle("Beta2")+
#   theme(axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_text(size=12,
#                                  color="black",
#                                  angle=0),
#         axis.text.x=element_text(size=12,
#                                  color="black",
#                                  angle=0),
#         plot.title = element_text(size=20,
#                                   hjust = 0.5))



ggplot(df_beta2_1,aes(x=beta2mcmc))+
  geom_histogram(aes(y=..density..,color="MCMC"),binwidth =0.005)+
  geom_line(data=df_beta2_2, aes(x=xSeqb1,y=newdensityBeta2,color="INFVB"),size=1.5)+
  geom_line(data=df_beta2_3, aes(x=xSeqb1,y=MFVB_beta2,color="MFVB"),size=1.5,linetype=5)+
  # geom_density(data = df_beta2_4, aes(x = beta2stan, fill = "STAN", color = "STAN"), alpha = 0.2, linetype = 1, size = 1.5) +
  scale_color_manual(name = "",values = c(MCMC = "black", INFVB = "blue",MFVB="red",STAN="yellow"))+
  coord_cartesian(xlim=c(0.85,1.1))+
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
foosigma2<-sample(x_sigma2,size = iter.thetaphi^2, replace = T, prob = ELBOvectorweights)
#MFVB
xSeqsigma2<-seq(0,1000, length.out=iter.thetaphi)


MFVB_sigma2<-dinvgamma(x=xSeqsigma2, shape = MFVB_iter[j,1] ,rate=MFVB_iter[j,2])
df_sigma2_1<-data.frame('sigma2mcmc'=parameterMatrix[-(1:burninperiod),3])
df_sigma2_2<-data.frame(foosigma2)
df_sigma2_3<-data.frame(MFVB_sigma2,xSeqsigma2)
# df_sigma2_4<-data.frame('sigma2stan'=parMat[,3])


ggplot(df_sigma2_1,aes(x=sigma2mcmc))+
  geom_histogram(aes(y=..density..,color="MCMC"),binwidth =20)+
  geom_density(data=df_sigma2_2, aes(x=foosigma2, color="INFVB"),size=1.5)+
  geom_line(data=df_sigma2_3, aes(x=xSeqsigma2,y=MFVB_sigma2,color="MFVB"),size=1.5,linetype=5)+
  # geom_density(data = df_sigma2_4, aes(x = sigma2stan, fill = "STAN", color = "STAN"), alpha = 0.2, linetype = 1, size = 1.5) +
  scale_color_manual(name = "",values = c(MCMC = "black", INFVB = "blue",MFVB="red",STAN="yellow"))+
  coord_cartesian(xlim=c(0,500))+
  # geom_vline(xintercept =sigma2, color="green",size=1.5)+
  theme(legend.position = "right")+
  ggtitle("sigma2")+
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
# Estimating the Eta  (Using MCMC, INFVB(phi), INFVB(phi,ssq), INLA)
################################################################################
################################################################################
# The real Data
eta_trainingdata_real <-X_trainingdata%*%beta+w_trainingdata
eta_testdata_real <-X_testdata%*%beta+w_testdata



#Choose the eta
ktheta=11
xSeq<-seq(-3,3,length.out=iter.thetaphi)
# INLA
ySeq_INLA<-dnorm(x = xSeq , mean = linear_predictor_mean[(0.8*n)+ktheta]  , sd = linear_predictor_sd[(0.8*n)+ktheta] )
# Sampling INFVBphi 
ySeq_VBfixssq_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[ktheta], sd=sqrt(cov_sample_eta_testdata_VB_all[ktheta,ktheta]))
# Sampling HMFVB 
ySeq_HMFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Pred[ktheta], sd=sqrt(cov_sample_eta_testdata_MFVB[ktheta,ktheta]))
#MCMC
Etamcmc<-eta_testdata_MCMC[ktheta,]
df_eta_MCMC<-data.frame(Etamcmc)
#Stan
# Etastan<-eta_testdata_MCMC_Stan[ktheta,]
# df_eta_STAN<-data.frame(Etamcmc)

df_eta_VBfixssq<-data.frame(ySeq_VBfixssq_sampling,xSeq)
df_eta_HMFVB<-data.frame(ySeq_HMFVB_sampling,xSeq)

df_eta_INLA<-data.frame(ySeq_INLA,xSeq)
colnames(df_eta_MCMC)[1]<-'Etamcmc'
# colnames(df_eta_STAN)[1]<-'Etastan'






ggplot() +
  geom_density(data = df_eta_MCMC, aes(x = Etamcmc, fill = "MCMC", color = "MCMC"), alpha = 0.6, linetype = "solid", size = 1.5) +
  geom_line(data = df_eta_VBfixssq, aes(x = xSeq, y = ySeq_VBfixssq_sampling, fill = "VBfixssq", color = "VBfixssq"), alpha = 0.6, linetype = "solid", size = 1) +
  geom_line(data = df_eta_HMFVB, aes(x = xSeq, y = ySeq_HMFVB_sampling, fill = "HMFVB", color = "HMFVB"), alpha = 0.6, linetype = "solid", size = 1) +
  geom_line(data = df_eta_INLA, aes(x = xSeq, y = ySeq_INLA, fill = "INLA", color = "INLA"), alpha = 0.6, linetype = "solid", size = 1) +
  # geom_density(data = df_eta_STAN, aes(x = Etastan, fill = "STAN", color = "STAN"), alpha = 0.2, linetype = 5, size = 1) +
  scale_fill_manual(
    name = "Methods",
    labels = c("MCMC ", "VB fixed ssq", "HMFVB", "STAN","INLA"),
    values = c(MCMC = "black", VBfixssq = "blue", HMFVB = "red", STAN = "yellow", INLA="orange")
  ) +
  scale_color_manual(
    name = "Methods",
    labels = c("MCMC ", "VB fixed ssq", "HMFVB", "STAN","INLA"),
    values = c(MCMC = "black", VBfixssq = "blue", HMFVB = "red", STAN = "yellow", INLA="orange")
  ) +
  geom_vline(xintercept = eta_testdata_real[ktheta], color = "green", size = 1.5, linetype = "solid") +
  # annotate("text", x = sigma2 + 0.1, y = 0.1, label = "σ²", color = "black", angle = 90, size = 4) +
  coord_cartesian(xlim = c(-1, 0)) +
  labs(x = "eta Values", y = "Density", title = "Posterior Distributions of Eta") +
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


#Choose the eta
ktheta=100
xSeq<-seq(-3,3,length.out=iter.thetaphi)
# INLA
ySeq_INLA<-dnorm(x = xSeq , mean = linear_predictor_mean[(0.8*n)+ktheta]  , sd = linear_predictor_sd[(0.8*n)+ktheta] )
# Sampling INFVBphi 
ySeq_VBfixssq_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[ktheta], sd=sqrt(cov_sample_eta_testdata_VB_all[ktheta,ktheta]))
# Sampling HMFVB 
ySeq_HMFVB_sampling <- dnorm(x=xSeq, mean=MFVB_mean_SA_Pred[ktheta], sd=sqrt(cov_sample_eta_testdata_MFVB[ktheta,ktheta]))
#MCMC
Etamcmc<-eta_testdata_MCMC[ktheta,]
df_eta_MCMC<-data.frame(Etamcmc)
#Stan
# Etastan<-eta_testdata_MCMC_Stan[ktheta,]
# df_eta_STAN<-data.frame(Etamcmc)

df_eta_VBfixssq<-data.frame(ySeq_VBfixssq_sampling,xSeq)
df_eta_HMFVB<-data.frame(ySeq_HMFVB_sampling,xSeq)

df_eta_INLA<-data.frame(ySeq_INLA,xSeq)
colnames(df_eta_MCMC)[1]<-'Etamcmc'
# colnames(df_eta_STAN)[1]<-'Etastan'






ggplot() +
  geom_density(data = df_eta_MCMC, aes(x = Etamcmc, fill = "MCMC", color = "MCMC"), alpha = 0.6, linetype = "solid", size = 1) +
  geom_line(data = df_eta_VBfixssq, aes(x = xSeq, y = ySeq_VBfixssq_sampling, fill = "VBfixssq", color = "VBfixssq"), alpha = 0.6, linetype = "solid", size = 1) +
  geom_line(data = df_eta_HMFVB, aes(x = xSeq, y = ySeq_HMFVB_sampling, fill = "HMFVB", color = "HMFVB"), alpha = 0.6, linetype = "solid", size = 1) +
  geom_line(data = df_eta_INLA, aes(x = xSeq, y = ySeq_INLA, fill = "INLA", color = "INLA"), alpha = 0.6, linetype = "solid", size = 1) +
  # geom_density(data = df_eta_STAN, aes(x = Etastan, fill = "STAN", color = "STAN"), alpha = 0.2, linetype = 5, size = 1) +
  scale_fill_manual(
    name = "Methods",
    labels = c("MCMC ", "VB fixed ssq", "HMFVB", "STAN","INLA"),
    values = c(MCMC = "black", VBfixssq = "blue", HMFVB = "red", STAN = "yellow", INLA="orange")
  ) +
  scale_color_manual(
    name = "Methods",
    labels = c("MCMC ", "VB fixed ssq", "HMFVB", "STAN","INLA"),
    values = c(MCMC = "black", VBfixssq = "blue", HMFVB = "red", STAN = "yellow", INLA="orange")
  ) +
  geom_vline(xintercept = eta_testdata_real[ktheta], color = "green", size = 1.5, linetype = "solid") +
  # annotate("text", x = sigma2 + 0.1, y = 0.1, label = "σ²", color = "black", angle = 90, size = 4) +
  coord_cartesian(xlim = c(-0.25, 0.5)) +
  labs(x = "eta Values", y = "Density", title = "Posterior Distributions of Eta") +
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



# dev.off()
################################################################################
################################################################################
################################################################################
################################################################################
#prob_testdata-Truth
#prob_testdata-Truth
df_prob_testdata<-data.frame(prob_testdata,locations_testdata)
head(df_prob_testdata)

ggplot(df_prob_testdata) +
  geom_point(aes(x = df_prob_testdata[, 2],
                 y = df_prob_testdata[, 3],
                 color = prob_testdata),
             shape = 7) +
  scale_color_gradientn(colours = rev(rainbow(3))) +  # <--- reversed color order
  theme(legend.position = "right",
        legend.title = element_blank()) +
  theme(axis.title.x = element_text(face = "plain",
                                    size = 15,
                                    hjust = 0.5,
                                    color = "black"),
        axis.title.y = element_text(face = "plain",
                                    size = 15,
                                    hjust = 0.5,
                                    color = "black"),
        axis.text.y = element_text(size = 12,
                                   color = "black",
                                   angle = 0),
        axis.text.x = element_text(size = 12,
                                   color = "black",
                                   angle = 0),
        plot.title = element_text(size = 20,
                                  hjust = 0.5))


################################################################################
# prob_testdata - MCMC
################################################################################
df_exp_prediction <- data.frame(exp_prediction, locations_testdata)
head(df_exp_prediction)

ggplot(df_exp_prediction) +
  geom_point(aes(x = df_exp_prediction[, 2],
                 y = df_exp_prediction[, 3],
                 color = exp_prediction_VB),
             shape = 7) +
  scale_color_gradientn(colours = rev(rainbow(3))) +  # reversed color scale
  theme(legend.position = "right",
        legend.title = element_blank()) +
  theme(axis.title.x = element_text(face = "plain", size = 15, hjust = 0.5, color = "black"),
        axis.title.y = element_text(face = "plain", size = 15, hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 12, color = "black", angle = 0),
        axis.text.x = element_text(size = 12, color = "black", angle = 0),
        plot.title = element_text(size = 20, hjust = 0.5))

################################################################################
# prob_testdata - INFVB
################################################################################
df_exp_prediction_VB <- data.frame(exp_prediction_VB, locations_testdata)
head(df_exp_prediction_VB)

ggplot(df_exp_prediction_VB) +
  geom_point(aes(x = df_exp_prediction_VB[, 2],
                 y = df_exp_prediction_VB[, 3],
                 color = exp_prediction_VB),
             shape = 7) +
  scale_color_gradientn(colours = rev(rainbow(3))) +  # reversed color scale
  theme(legend.position = "right",
        legend.title = element_blank()) +
  theme(axis.title.x = element_text(face = "plain", size = 15, hjust = 0.5, color = "black"),
        axis.title.y = element_text(face = "plain", size = 15, hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 12, color = "black", angle = 0),
        axis.text.x = element_text(size = 12, color = "black", angle = 0),
        plot.title = element_text(size = 20, hjust = 0.5))

################################################################################
# prob_testdata - MFVB
################################################################################
df_exp_prediction_MFVB <- data.frame(exp_prediction_MFVB, locations_testdata)
head(df_exp_prediction_MFVB)

ggplot(df_exp_prediction_MFVB) +
  geom_point(aes(x = df_exp_prediction_MFVB[, 2],
                 y = df_exp_prediction_MFVB[, 3],
                 color = exp_prediction_VB),
             shape = 7) +
  scale_color_gradientn(colours = rev(rainbow(3))) +  # reversed color scale
  theme(legend.position = "right",
        legend.title = element_blank()) +
  theme(axis.title.x = element_text(face = "plain", size = 15, hjust = 0.5, color = "black"),
        axis.title.y = element_text(face = "plain", size = 15, hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 12, color = "black", angle = 0),
        axis.text.x = element_text(size = 12, color = "black", angle = 0),
        plot.title = element_text(size = 20, hjust = 0.5))




################################################################################
################################################################################
## RMSPE FINAL RESULT Comparison between MCMC vs VB 
################################################################################
################################################################################
# obs_RMSPE # RMSPE - MCMC
# obs_RMSPE_MFVB # RMSPE - MFVB
# obs_RMSPE_VB # RMSPE - INFVB

aucVal
aucVal_MFVB
aucVal_VB

ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)

