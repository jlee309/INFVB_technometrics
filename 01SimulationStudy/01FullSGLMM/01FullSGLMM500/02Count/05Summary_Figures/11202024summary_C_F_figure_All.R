################################################################################
################################################################################
# 11202024 Figures
################################################################################
################################################################################

setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/02Count/05Summary_Figures")
getwd()

dev.off()
library(mvtnorm)
library(cli);library(viridis);library(fields);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm)
library(invgamma);library(pROC);library(Matrix);library(emulator);library(ggplot2)
getwd()
################################################################################
################################################################################
# rm(list=ls())
# load(file="../01DataGeneration/C_phi1_Gen500_1.RData")
# load(file="../02MCMC/C_phi1_Gen500_1_MCMC.RData")
# load(file="../04INFVBfixssqphi_ParallelComputing/C_phi1_Gen500_1_VBphissqphiParallel.RData") 
# load(file="../03INFVBfixphi_ParallelComputing/C_phi1_Gen500_1_VBphiParallel.RData")
# load(file="../06INLA/C_phi1_Gen500_1_INLA.RData")
# load(file="../10Stan/C_phi1_Gen500_1_MCMC_Stan.RData")
# load(file="../08INFVBfixphi_SA/SA_C_phi1_Gen500_1_Full_INFVBphi.RData")
# load(file="../09INFVBfixssqphi_SA/SA_C_phi1_Gen500_1_Full_INFVBssqphi.RData")

# rm(list=ls())
# load(file="../01DataGeneration/C_phi2_Gen500_1.RData")
# load(file="../02MCMC/C_phi2_Gen500_1_MCMC.RData")
# load(file="../04INFVBfixssqphi_ParallelComputing/C_phi2_Gen500_1_VBphissqphiParallel.RData")
# load(file="../03INFVBfixphi_ParallelComputing/C_phi2_Gen500_1_VBphiParallel.RData")
# load(file="../06INLA/C_phi2_Gen500_1_INLA.RData")
# load(file="../10Stan/C_phi2_Gen500_1_MCMC_Stan.RData")
# load(file="../08INFVBfixphi_SA/SA_C_phi2_Gen500_1_Full_INFVBphi.RData")
# load(file="../09INFVBfixssqphi_SA/SA_C_phi2_Gen500_1_Full_INFVBssqphi.RData")
# 
# rm(list=ls())
# load(file="../01DataGeneration/C_phi3_Gen500_1.RData")
# load(file="../02MCMC/C_phi3_Gen500_1_MCMC.RData")
# load(file="../04INFVBfixssqphi_ParallelComputing/C_phi3_Gen500_1_VBphissqphiParallel.RData") 
# load(file="../03INFVBfixphi_ParallelComputing/C_phi3_Gen500_1_VBphiParallel.RData")
# load(file="../06INLA/C_phi3_Gen500_1_INLA.RData")
# load(file="../10Stan/C_phi3_Gen500_1_MCMC_Stan.RData")
# load(file="../08INFVBfixphi_SA/SA_C_phi3_Gen500_1_Full_INFVBphi.RData")
# load(file="../09INFVBfixssqphi_SA/SA_C_phi3_Gen500_1_Full_INFVBssqphi.RData")
# 
# rm(list=ls())
load(file="../01DataGeneration/C_phi4_Gen500_1.RData")
load(file="../02MCMC/C_phi4_Gen500_1_MCMC.RData")
load(file="../04INFVBfixssqphi_ParallelComputing/C_phi4_Gen500_1_VBphissqphiParallel.RData")
load(file="../03INFVBfixphi_ParallelComputing/C_phi4_Gen500_1_VBphiParallel.RData")
load(file="../06INLA/C_phi4_Gen500_1_INLA.RData")
# load(file="../10Stan/C_phi4_Gen500_1_MCMC_Stan.RData")
load(file="../08INFVBfixphi_SA/SA_C_phi4_Gen500_1_Full_INFVBphi.RData")
load(file="../09INFVBfixssqphi_SA/SA_C_phi4_Gen500_1_Full_INFVBssqphi.RData")


################################################################################
################################################################################
# Estimating the parameter (beta,ssq,phi) (Using MCMC, INFVB(phi), INFVB(phi,ssq), Stan)
################################################################################
################################################################################
################################################################################
# Beta 1
################################################################################
dev.off()
xSeqb1<-seq(-2,2, length.out=length(newdensityBeta1))
xSeqb1ssqphi<-seq(-3,3, length.out=length(newdensityBeta1ssqphi))

beta1mcmc<-parameterMatrix[-(1:burninperiod),1]
# beta1stan<-parMat[,1]
df_beta1_1<-data.frame(beta1mcmc)
df_beta1_2<-data.frame(newdensityBeta1,xSeqb1)
df_beta1_3<-data.frame(newdensityBeta1ssqphi,xSeqb1ssqphi)
# df_beta1_4<-data.frame(beta1stan)
colnames(df_beta1_1)[1]<-'MCMC_beta1'
# colnames(df_beta1_4)[1]<-'Stan_beta1'

ggplot() +
  geom_density(data = df_beta1_1, aes(x = MCMC_beta1, fill = "MCMC", color = "MCMC"), alpha = 0.2, linetype = "solid", size = 1.5) +
  geom_line(data = df_beta1_2, aes(x = xSeqb1, y = newdensityBeta1, fill = "VBfixphi", color = "VBfixphi"), alpha = 0.6, linetype = "solid", size = 1.5) +
  geom_line(data = df_beta1_3, aes(x = xSeqb1ssqphi, y = newdensityBeta1ssqphi, fill = "VBfixssqphi", color = "VBfixssqphi"), alpha = 0.6, linetype = 5, size = 1.5) +
  # geom_density(data = df_beta1_4, aes(x = Stan_beta1, fill = "STAN", color = "STAN"), alpha = 0.2, linetype = "solid", size = 1.5) +
  scale_fill_manual(
    name = "Methods",
    labels = c("MCMC Samples", "VB fixed φ", "VB fixed σ²φ", "STAN Samples"),
    values = c(MCMC = "darkorange", VBfixphi = "purple", VBfixssqphi = "cyan", STAN = "limegreen")
  ) +
  scale_color_manual(
    name = "Methods",
    labels = c("MCMC Samples", "VB fixed φ", "VB fixed σ²φ", "STAN Samples"),
    values = c(MCMC = "darkorange", VBfixphi = "purple", VBfixssqphi = "cyan", STAN = "limegreen")
  ) +
  geom_vline(xintercept = beta[1], color = "black", size = 1.5, linetype = "dashed") +
  # annotate("text", x = sigma2 + 0.1, y = 0.1, label = "σ²", color = "black", angle = 90, size = 4) +
  coord_cartesian(xlim = c(0.7, 1.3)) +
  labs(x = "Beta1 Values", y = "Density", title = "Posterior Distributions of Beta1") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )


################################################################################
# Beta 2
################################################################################
xSeqb2<-seq(-2,2, length.out=length(newdensityBeta2))
xSeqb2ssqphi<-seq(-3,3, length.out=length(newdensityBeta2ssqphi))

beta2mcmc<-parameterMatrix[-(1:burninperiod),2]
# beta2stan<-parMat[,2]
df_beta2_1<-data.frame(beta2mcmc)
df_beta2_2<-data.frame(newdensityBeta2,xSeqb2)
df_beta2_3<-data.frame(newdensityBeta2ssqphi,xSeqb2ssqphi)
# df_beta2_4<-data.frame(beta2stan)
colnames(df_beta2_1)[1]<-'MCMC_beta2'
# colnames(df_beta2_4)[1]<-'Stan_beta2'

ggplot() +
  geom_density(data = df_beta2_1, aes(x = MCMC_beta2, fill = "MCMC", color = "MCMC"), alpha = 0.2, linetype = "solid", size = 1.5) +
  geom_line(data = df_beta2_2, aes(x = xSeqb2, y = newdensityBeta2, fill = "VBfixphi", color = "VBfixphi"), alpha = 0.6, linetype = "solid", size = 1.5) +
  geom_line(data = df_beta2_3, aes(x = xSeqb2ssqphi, y = newdensityBeta2ssqphi, fill = "VBfixssqphi", color = "VBfixssqphi"), alpha = 0.6, linetype = 5, size = 1.5) +
  # geom_density(data = df_beta2_4, aes(x = Stan_beta2, fill = "STAN", color = "STAN"), alpha = 0.2, linetype = "solid", size = 1.5) +
  scale_fill_manual(
    name = "Methods",
    labels = c("MCMC Samples", "VB fixed φ", "VB fixed σ²φ", "STAN Samples"),
    values = c(MCMC = "darkorange", VBfixphi = "purple", VBfixssqphi = "cyan", STAN = "limegreen")
  ) +
  scale_color_manual(
    name = "Methods",
    labels = c("MCMC Samples", "VB fixed φ", "VB fixed σ²φ", "STAN Samples"),
    values = c(MCMC = "darkorange", VBfixphi = "purple", VBfixssqphi = "cyan", STAN = "limegreen")
  ) +
  geom_vline(xintercept = beta[2], color = "black", size = 1.5, linetype = "dashed") +
  # annotate("text", x = sigma2 + 0.1, y = 0.1, label = "σ²", color = "black", angle = 90, size = 4) +
  coord_cartesian(xlim = c(0.75, 1.5)) +
  labs(x = "Beta2 Values", y = "Density", title = "Posterior Distributions of Beta2") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )





################################################################################
## Sigma2
################################################################################

set.seed(123)
foosigma2<-sample(x_sigma2,size = iter.thetaphissqphi, replace = T, prob = ELBOsigma2)
foosigma2<-data.frame(foosigma2)

xSeqsigma2<-seq(0.1,6, length.out=length(newdensitySigma2))
foosigma2fixphi<-data.frame(xSeqsigma2,newdensitySigma2)

MCMC_sigma2<-parameterMatrix[-(1:burninperiod),4]
MCMC_sigma2<-data.frame(MCMC_sigma2)

# Stan_sigma2<-parMat[,4]
# Stan_sigma2<-data.frame(Stan_sigma2)

colnames(foosigma2)[1]<-'VBfixssqphi'
colnames(foosigma2fixphi)[1]<-'xSeqsigma2';colnames(foosigma2fixphi)[2]<-'newdensitySigma2'
colnames(MCMC_sigma2)[1]<-'MCMC_sigma2'
# colnames(Stan_sigma2)[1]<-'Stan_sigma2'

ggplot()+
  geom_density(data=MCMC_sigma2, aes(x=MCMC_sigma2, fill="MCMC",color="MCMC"),alpha=0.2, linetype="solid",size=1.5)+
  geom_line(data=foosigma2fixphi, aes(x=xSeqsigma2,y=newdensitySigma2,fill="VBfixphi",color="VBfixphi"),alpha=0.6,linetype="solid",size=1.5)+
  geom_density(data=foosigma2, aes(x=VBfixssqphi,fill="VBfixssqphi",color="VBfixssqphi"),alpha=0.4,linetype="solid",size=1.5)+
  # geom_density(data=Stan_sigma2, aes(x=Stan_sigma2, fill="STAN",color="STAN"),alpha=0.2, linetype="solid",size=1.5)+
  scale_fill_manual(
    name = "Methods",
    labels = c("MCMC Samples", "VB fixed φ", "VB fixed σ²φ", "STAN Samples"),
    values = c(MCMC = "darkorange", VBfixphi = "purple", VBfixssqphi = "cyan", STAN = "limegreen")
  ) +
  scale_color_manual(
    name = "Methods",
    labels = c("MCMC Samples", "VB fixed φ", "VB fixed σ²φ", "STAN Samples"),
    values = c(MCMC = "darkorange", VBfixphi = "purple", VBfixssqphi = "cyan", STAN = "limegreen")
  ) +
  geom_vline(xintercept = sigma2, color = "black", size = 1.5, linetype = "dashed") +
  # annotate("text", x = sigma2 + 0.1, y = 0.1, label = "σ²", color = "black", angle = 90, size = 4) +
  coord_cartesian(xlim = c(0, 2.5)) +
  labs(x = "Sigma2 Values", y = "Density", title = "Posterior Distributions of sigma2") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )


################################################################################
## Phi
################################################################################
set.seed(1234)
foo<-sample(thetaphi,size = iter.thetaphi, replace = T, prob = ELBOvectorweights )
foo<-data.frame(foo) #VBfixphi

foophi<-sample(thetaphissqphi,size = iter.thetaphissqphi, replace = T, prob = ELBOphi )
foophi<-data.frame(foophi) # VBfixssqphi

MCMC_phi<-parameterMatrix[-(1:burninperiod),3] #MCMCphi
MCMC_phi<-data.frame(MCMC_phi)

# Stan_phi<-parMat[,3] #MCMCphi
# Stan_phi<-data.frame(Stan_phi)


colnames(foo)[1]<-'VBfixphi'
colnames(foophi)[1]<-'VBfixssqphi'
colnames(MCMC_phi)[1]<-'MCMCphi'
# colnames(Stan_phi)[1]<-'Stanphi'

ggplot()+
  geom_density(data=MCMC_phi, aes(x=MCMCphi, fill="MCMC",color="MCMC"),alpha=0.2, linetype="solid",size=1.5)+
  geom_density(data=foo, aes(x=VBfixphi,fill="VBfixphi",color="VBfixphi"),alpha=0.3,linetype="solid",size=1.5)+
  geom_density(data=foophi, aes(x=VBfixssqphi, fill="VBfixssqphi",color="VBfixssqphi"),alpha=0.4,linetype="solid",size=1.5)+
  # geom_density(data=Stan_phi, aes(x=Stanphi, fill="STAN",color="STAN"),alpha=0.4,linetype="solid",size=1.5)+
  scale_fill_manual(
    name = "Methods",
    labels = c("MCMC Samples", "VB fixed φ", "VB fixed σ²φ", "STAN Samples"),
    values = c(MCMC = "darkorange", VBfixphi = "purple", VBfixssqphi = "cyan", STAN = "limegreen")
  ) +
  scale_color_manual(
    name = "Methods",
    labels = c("MCMC Samples", "VB fixed φ", "VB fixed σ²φ", "STAN Samples"),
    values = c(MCMC = "darkorange", VBfixphi = "purple", VBfixssqphi = "cyan", STAN = "limegreen")
  ) +
  geom_vline(xintercept = phi, color = "black", size = 1.5, linetype = "dashed") +
  # annotate("text", x = sigma2 + 0.1, y = 0.1, label = "σ²", color = "black", angle = 90, size = 4) +
  coord_cartesian(xlim = c(0, sqrt(2))) +
  labs(x = "phi Values", y = "Density", title = "Posterior Distributions of phi") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

################################################################################
################################################################################
# Estimating the Eta  (Using MCMC, INFVB(phi), INFVB(phi,ssq), INLA)
################################################################################
################################################################################
# The real Data
eta_trainingdata_real <-X_trainingdata%*%beta+w_trainingdata
eta_testdata_real <-X_testdata%*%beta+w_testdata

#Choose the eta
ktheta=27

xSeq<-seq(-3,3,length.out=iter.thetaphi)
# INLA
ySeq_INLA<-dnorm(x = xSeq , mean = field_pred_mean[ktheta] , sd = field_pred_sd[ktheta])
# Sampling INFVBphi 
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[ktheta], sd=sqrt(cov_sample_eta_testdata_VB_all[ktheta,ktheta]))
# Sampling INFVBssqphi 
ySeq_VBfixssqphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VBssqphi[ktheta], sd=sqrt(cov_sample_eta_testdata_VBssqphi_all[ktheta,ktheta]))

Etamcmc<-eta_testdata_MCMC[ktheta,-(1:burninperiod)]
df_eta_MCMC<-data.frame(Etamcmc)

# Etastan<-eta_testdata_MCMC_Stan[ktheta,]
# df_eta_STAN<-data.frame(Etamcmc)

df_eta_VBfixphi<-data.frame(ySeq_VBfixphi_sampling,xSeq)
df_eta_VBfixssqphi<-data.frame(ySeq_VBfixssqphi_sampling,xSeq)

df_eta_INLA<-data.frame(ySeq_INLA,xSeq)
colnames(df_eta_MCMC)[1]<-'Etamcmc'
# colnames(df_eta_STAN)[1]<-'Etastan'


ggplot() +
  geom_density(data = df_eta_MCMC, aes(x = Etamcmc, fill = "MCMC", color = "MCMC"), alpha = 0.2, linetype = "solid", size = 1.5) +
  geom_line(data = df_eta_VBfixphi, aes(x = xSeq, y = ySeq_VBfixphi_sampling, fill = "VBfixphi", color = "VBfixphi"), alpha = 0.6, linetype = "solid", size = 1.5) +
  geom_line(data = df_eta_VBfixssqphi, aes(x = xSeq, y = ySeq_VBfixssqphi_sampling, fill = "VBfixssqphi", color = "VBfixssqphi"), alpha = 0.6, linetype = 5, size = 1.5) +
  geom_line(data = df_eta_INLA, aes(x = xSeq, y = ySeq_INLA, fill = "INLA", color = "INLA"), alpha = 0.6, linetype = "solid", size = 1.5) +
  # geom_density(data = df_eta_STAN, aes(x = Etastan, fill = "STAN", color = "STAN"), alpha = 0.2, linetype = 5, size = 1.5) +
  scale_fill_manual(
    name = "Methods",
    labels = c("MCMC ", "VB fixed φ", "VB fixed σ²φ", "STAN","INLA"),
    values = c(MCMC = "darkorange", VBfixphi = "purple", VBfixssqphi = "cyan", STAN = "limegreen", INLA="red")
  ) +
  scale_color_manual(
    name = "Methods",
    labels = c("MCMC ", "VB fixed φ", "VB fixed σ²φ", "STAN","INLA"),
    values = c(MCMC = "darkorange", VBfixphi = "purple", VBfixssqphi = "cyan", STAN = "limegreen",INLA="red")
  ) +
  geom_vline(xintercept = eta_testdata_real[ktheta], color = "black", size = 1.5, linetype = "dashed") +
  # annotate("text", x = sigma2 + 0.1, y = 0.1, label = "σ²", color = "black", angle = 90, size = 4) +
  coord_cartesian(xlim = c(-2.5, 0.3)) +
  labs(x = "eta Values", y = "Density", title = "Posterior Distributions of Eta") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )


################################################################################
################################################################################
################################################################################
#Choose the eta
ktheta=14
xSeq<-seq(-3,3,length.out=iter.thetaphi)
# INLA
ySeq_INLA<-dnorm(x = xSeq , mean = field_pred_mean[ktheta] , sd = field_pred_sd[ktheta])
# Sampling INFVBphi 
ySeq_VBfixphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VB[ktheta], sd=sqrt(cov_sample_eta_testdata_VB_all[ktheta,ktheta]))
# Sampling INFVBssqphi 
ySeq_VBfixssqphi_sampling <- dnorm(x=xSeq, mean=mu_sample_pred_VBssqphi[ktheta], sd=sqrt(cov_sample_eta_testdata_VBssqphi_all[ktheta,ktheta]))

Etamcmc<-eta_testdata_MCMC[ktheta,-(1:burninperiod)]
df_eta_MCMC<-data.frame(Etamcmc)

# Etastan<-eta_testdata_MCMC_Stan[ktheta,]
# df_eta_STAN<-data.frame(Etamcmc)

df_eta_VBfixphi<-data.frame(ySeq_VBfixphi_sampling,xSeq)
df_eta_VBfixssqphi<-data.frame(ySeq_VBfixssqphi_sampling,xSeq)

df_eta_INLA<-data.frame(ySeq_INLA,xSeq)
colnames(df_eta_MCMC)[1]<-'Etamcmc'
# colnames(df_eta_STAN)[1]<-'Etastan'


ggplot() +
  geom_density(data = df_eta_MCMC, aes(x = Etamcmc, fill = "MCMC", color = "MCMC"), alpha = 0.2, linetype = "solid", size = 1.5) +
  geom_line(data = df_eta_VBfixphi, aes(x = xSeq, y = ySeq_VBfixphi_sampling, fill = "VBfixphi", color = "VBfixphi"), alpha = 0.6, linetype = "solid", size = 1.5) +
  geom_line(data = df_eta_VBfixssqphi, aes(x = xSeq, y = ySeq_VBfixssqphi_sampling, fill = "VBfixssqphi", color = "VBfixssqphi"), alpha = 0.6, linetype = 5, size = 1.5) +
  geom_line(data = df_eta_INLA, aes(x = xSeq, y = ySeq_INLA, fill = "INLA", color = "INLA"), alpha = 0.6, linetype = "solid", size = 1.5) +
  # geom_density(data = df_eta_STAN, aes(x = Etastan, fill = "STAN", color = "STAN"), alpha = 0.2, linetype = 5, size = 1.5) +
  scale_fill_manual(
    name = "Methods",
    labels = c("MCMC ", "VB fixed φ", "VB fixed σ²φ", "STAN","INLA"),
    values = c(MCMC = "darkorange", VBfixphi = "purple", VBfixssqphi = "cyan", STAN = "limegreen", INLA="red")
  ) +
  scale_color_manual(
    name = "Methods",
    labels = c("MCMC ", "VB fixed φ", "VB fixed σ²φ", "STAN","INLA"),
    values = c(MCMC = "darkorange", VBfixphi = "purple", VBfixssqphi = "cyan", STAN = "limegreen",INLA="red")
  ) +
  geom_vline(xintercept = eta_testdata_real[ktheta], color = "black", size = 1.5, linetype = "dashed") +
  # annotate("text", x = sigma2 + 0.1, y = 0.1, label = "σ²", color = "black", angle = 90, size = 4) +
  coord_cartesian(xlim = c(-1.0, 1.5)) +
  labs(x = "eta Values", y = "Density", title = "Posterior Distributions of Eta") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
