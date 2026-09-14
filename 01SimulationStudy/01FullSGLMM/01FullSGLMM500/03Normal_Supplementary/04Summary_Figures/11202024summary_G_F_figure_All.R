################################################################################
################################################################################
# SUMMARY RESULTS FOR NORMAL when n=500 Supplementary Figure 2 phi=0.4
################################################################################
################################################################################
rm(list=ls())
setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/03Normal_Supplementary/04Summary_Figures")
getwd()
dev.off()

# load(file="../01DataGeneration/G_phi2_Gen500_1.RData")
# load(file="../02MCMC/G_phi2_Gen500_1_MCMC.RData")
# load(file="../03INFVBfixphi_ParallelComputing/G_phi2_Gen500_1_VBphiParallel.RData")
# load(file="../05INFVBfixphi_SA/SA_G_phi2_Gen500_1_Full_INFVBphi.RData")
# load(file="../10Stan/G_phi2_Gen500_1_MCMC_Stan.RData")

load(file="../01DataGeneration/G_phi4_Gen500_1.RData")
load(file="../02MCMC/G_phi4_Gen500_1_MCMC.RData")
load(file="../03INFVBfixphi_ParallelComputing/G_phi4_Gen500_1_VBphiParallel.RData")
load(file="../05INFVBfixphi_SA/SA_G_phi4_Gen500_1_Full_INFVBphi.RData")
load(file="../10Stan/G_phi4_Gen500_1_MCMC_Stan.RData")

library(cli);library(viridis);library(fields);library(pgdraw);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm);
library(invgamma);library(pROC);library(Matrix);library(emulator);
library(ggplot2);library(ggthemes)



################################################################################
################################################################################
################################################################################
################################################################################
## ggplot
################################################################################
################################################################################
################################################################################
################################################################################
## Beta 1
################################################################################
library(ggplot2)
beta1mcmc<-parameterMatrix[-(1:burninperiod),1]
df_beta1_1<-data.frame(beta1mcmc)
df_beta1_2<-data.frame(newdensityBeta1,xSeqb1)
colnames(df_beta1_1)[1]<-'MCMC_beta1'

ggplot()+
  geom_density(data=df_beta1_1, aes(x=MCMC_beta1, fill="MCMC",color="MCMC"),alpha=0.2, linetype="solid",size=1.5)+
  geom_line(data=df_beta1_2, aes(x=xSeqb1,y=newdensityBeta1,fill="VBfixphi",color="VBfixphi"),alpha=0.6,linetype="solid",size=1.5)+
scale_fill_manual(
  name = "Methods",
  labels = c("MCMC Samples", "VB fixed φ"),
  values = c(MCMC = "darkorange", VBfixphi = "purple")
) +
  scale_color_manual(
    name = "Methods",
    labels = c("MCMC Samples", "VB fixed φ"),
    values = c(MCMC = "darkorange", VBfixphi = "purple")
  ) +
  geom_vline(xintercept = beta[1], color = "black", size = 1.5, linetype = "dashed") +
  # annotate("text", x = sigma2 + 0.1, y = 0.1, label = "σ²", color = "black", angle = 90, size = 4) +
  coord_cartesian(xlim = c(0.9, 1.1)) +
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
## Beta 2
################################################################################
beta2mcmc<-parameterMatrix[-(1:burninperiod),2]
df_beta2_1<-data.frame(beta2mcmc)
df_beta2_2<-data.frame(newdensityBeta2,xSeqb2)
colnames(df_beta2_1)[1]<-'MCMC_beta2'


ggplot()+
  geom_density(data=df_beta2_1, aes(x=MCMC_beta2, fill="MCMC",color="MCMC"),alpha=0.2, linetype="solid",size=1.5)+
  geom_line(data=df_beta2_2, aes(x=xSeqb1,y=newdensityBeta2,fill="VBfixphi",color="VBfixphi"),alpha=0.6,linetype="solid",size=1.5)+
scale_fill_manual(
  name = "Methods",
  labels = c("MCMC Samples", "VB fixed φ"),
  values = c(MCMC = "darkorange", VBfixphi = "purple")
) +
  scale_color_manual(
    name = "Methods",
    labels = c("MCMC Samples", "VB fixed φ"),
    values = c(MCMC = "darkorange", VBfixphi = "purple")
  ) +
  geom_vline(xintercept = beta[2], color = "black", size = 1.5, linetype = "dashed") +
  # annotate("text", x = sigma2 + 0.1, y = 0.1, label = "σ²", color = "black", angle = 90, size = 4) +
  coord_cartesian(xlim = c(0.95, 1.1)) +
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
foosigma2fixphi<-data.frame(xSeqsigma2,newdensitySigma2)
MCMC_sigma2<-parameterMatrix[-(1:burninperiod),3]
MCMC_sigma2<-data.frame(MCMC_sigma2)
colnames(foosigma2fixphi)[1]<-'xSeqsigma2';colnames(foosigma2fixphi)[2]<-'newdensitySigma2'
colnames(MCMC_sigma2)[1]<-'MCMC_sigma2'
head(foosigma2fixphi)
head(MCMC_sigma2)

# ggplot()+
#   geom_density(data=MCMC_sigma2, aes(x=MCMC_sigma2, fill="MCMC",color="MCMC"),alpha=0.2, linetype="solid",size=1.5)+
#   geom_line(data=foosigma2fixphi, aes(x=xSeqsigma2,y=newdensitySigma2,fill="VBfixphi",color="VBfixphi"),alpha=0.6,linetype="solid",size=1.5)+
#   scale_fill_manual(name = "", values = c(MCMC = "red", VBfixphi = "blue"))+
#   scale_color_manual(name = "",values = c(MCMC = "red", VBfixphi = "blue"))+
#   geom_vline(xintercept =sigma2, color="black",size=1)+
#   coord_cartesian(xlim=c(0,4))+
#   theme(legend.position = "right")+
#   ggtitle("Sigma2")+
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
# 


ggplot()+
  geom_density(data=MCMC_sigma2, aes(x=MCMC_sigma2, fill="MCMC",color="MCMC"),alpha=0.2, linetype="solid",size=1.5)+
  geom_line(data=foosigma2fixphi, aes(x=xSeqsigma2,y=newdensitySigma2,fill="VBfixphi",color="VBfixphi"),alpha=0.6,linetype="solid",size=1.5)+
  scale_fill_manual(
    name = "Methods",
    labels = c("MCMC Samples", "VB fixed φ"),
    values = c(MCMC = "darkorange", VBfixphi = "purple")
  ) +
  scale_color_manual(
    name = "Methods",
    labels = c("MCMC Samples", "VB fixed φ"),
    values = c(MCMC = "darkorange", VBfixphi = "purple")
  ) +
  geom_vline(xintercept = sigma2, color = "black", size = 1.5, linetype = "dashed") +
  # annotate("text", x = sigma2 + 0.1, y = 0.1, label = "σ²", color = "black", angle = 90, size = 4) +
  coord_cartesian(xlim = c(0, 2)) +
  labs(x = "Sigma2 Values", y = "Density", title = "Posterior Distributions of Sigma2") +
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
set.seed(20252)
foo<-sample(thetaphi,size = iter.thetaphi, replace = T, prob = ELBOvectorweights )
foo<-data.frame(foo) #VBfixphi
MCMC_phi<-parameterMatrix[-(1:burninperiod),4] #MCMCphi
MCMC_phi<-data.frame(MCMC_phi)

colnames(foo)[1]<-'VBfixphi'
colnames(MCMC_phi)[1]<-'MCMCphi'

# ggplot()+
#   geom_density(data=MCMC_phi, aes(x=MCMCphi, fill="MCMC",color="MCMC"),alpha=0.3, linetype="solid",size=1.5)+
#   geom_density(data=foo, aes(x=VBfixphi,fill="VBfixphi",color="VBfixphi"),alpha=0.05,linetype="solid",size=1.5)+
#   scale_fill_manual(name = "", values = c(MCMC = "red", VBfixphi = "blue"))+
#   scale_color_manual(name = "",values = c(MCMC = "red", VBfixphi = "blue")) +
#   geom_vline(xintercept =phi, color="black",size=2)+
#   coord_cartesian(xlim=c(0,1))+
#   theme(legend.position = "right")+
#   ggtitle("Phi")+
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
# 



ggplot()+
  geom_density(data=MCMC_phi, aes(x=MCMCphi, fill="MCMC",color="MCMC"),alpha=0.3, linetype="solid",size=1.5)+
  geom_density(data=foo, aes(x=VBfixphi,fill="VBfixphi",color="VBfixphi"),alpha=0.05,linetype="solid",size=1.5)+
scale_fill_manual(
  name = "Methods",
  labels = c("MCMC Samples", "VB fixed φ"),
  values = c(MCMC = "darkorange", VBfixphi = "purple")
) +
  scale_color_manual(
    name = "Methods",
    labels = c("MCMC Samples", "VB fixed φ"),
    values = c(MCMC = "darkorange", VBfixphi = "purple")
  ) +
  geom_vline(xintercept = phi, color = "black", size = 1.5, linetype = "dashed") +
  # annotate("text", x = sigma2 + 0.1, y = 0.1, label = "σ²", color = "black", angle = 90, size = 4) +
  coord_cartesian(xlim = c(0, 1)) +
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
## RMSPE FINAL RESULT Comparison between MCMC vs VB 
################################################################################
################################################################################
Z_RMSPE_MCMC # RMSPE - MCMC
Z_RMSPE_VB # RMSPE - VB fixphi
ptFinal # Walltime - MCMC
VB_ptFinal # Walltime - VB fixphi
ptFinal/VB_ptFinal # Computational Speedup 

