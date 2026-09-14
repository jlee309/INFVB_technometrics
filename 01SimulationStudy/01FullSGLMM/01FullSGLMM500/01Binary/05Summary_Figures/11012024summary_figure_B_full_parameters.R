
setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/01Binary/05Summary_Figures")
getwd()

rm(list=ls())
dev.off()
library(mvtnorm)
library(cli);library(viridis);library(fields);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm)
library(invgamma);library(pROC);library(Matrix);library(emulator);library(ggplot2)


load(file="../01DataGeneration/B_phi1_Gen500_1.RData")
load(file="../02MCMC/B_phi1_Gen500_1_MCMC.RData")
load(file="../04INFVBfixssqphi_ParallelComputing/B_phi1_Gen500_1_VBphissqphiParallel.RData") #elbovectorweight has some problem here
load(file="../03INFVBfixphi_ParallelComputing/B_phi1_Gen500_1_VBphiParallel.RData")
load(file="../06INLA/B_phi1_Gen500_1_INLA.RData")
load(file="../10Stan/B_phi1_Gen500_1_MCMC_Stan.RData")


################################################################################
################################################################################




################################################################################
################################################################################
################################################################################
################################################################################
# ggplot
################################################################################
################################################################################
################################################################################
################################################################################
# Beta 1
################################################################################
xSeqb1<-seq(-3,3, length.out=length(newdensityBeta1))
xSeqb1ssqphi<-seq(-3,3, length.out=length(newdensityBeta1ssqphi))

beta1mcmc<-parameterMatrix[-(1:burninperiod),1]
df_beta1_1<-data.frame(beta1mcmc)
df_beta1_2<-data.frame(newdensityBeta1,xSeqb1)
df_beta1_3<-data.frame(newdensityBeta1ssqphi,xSeqb1ssqphi)
colnames(df_beta1_1)[1]<-'MCMC_beta1'

ggplot()+
  geom_density(data=df_beta1_1, aes(x=MCMC_beta1, fill="MCMC",color="MCMC"),alpha=0.2, linetype="solid",size=1.5)+
  geom_line(data=df_beta1_2, aes(x=xSeqb1,y=newdensityBeta1,fill="VBfixphi",color="VBfixphi"),alpha=0.6,linetype="solid",size=1.5)+
  geom_line(data=df_beta1_3, aes(x=xSeqb1ssqphi,y=newdensityBeta1ssqphi,fill="VBfixssqphi",color="VBfixssqphi"),alpha=0.6,linetype="solid",size=1.5)+
  scale_fill_manual(name = "", values = c(MCMC = "red", VBfixphi = "blue",VBfixssqphi="green"))+
  scale_color_manual(name = "",values = c(MCMC = "red", VBfixphi = "blue",VBfixssqphi="green"))+
  geom_vline(xintercept =sigma2, color="black",size=1.5)+
  coord_cartesian(xlim=c(0,2))+
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
xSeqb2<-seq(-3,3, length.out=length(newdensityBeta2))
xSeqb2ssqphi<-seq(-3,3, length.out=length(newdensityBeta2ssqphi))

beta2mcmc<-parameterMatrix[-(1:burninperiod),2]
df_beta2_1<-data.frame(beta2mcmc)
df_beta2_2<-data.frame(newdensityBeta2,xSeqb2)
df_beta2_3<-data.frame(newdensityBeta2ssqphi,xSeqb2ssqphi)
colnames(df_beta2_1)[1]<-'MCMC_beta2'

ggplot()+
  geom_density(data=df_beta2_1, aes(x=MCMC_beta2, fill="MCMC",color="MCMC"),alpha=0.2, linetype="solid",size=1.5)+
  geom_line(data=df_beta2_2, aes(x=xSeqb1,y=newdensityBeta2,fill="VBfixphi",color="VBfixphi"),alpha=0.6,linetype="solid",size=1.5)+
  geom_line(data=df_beta2_3, aes(x=xSeqb1ssqphi,y=newdensityBeta2ssqphi,fill="VBfixssqphi",color="VBfixssqphi"),alpha=0.6,linetype="solid",size=1.5)+
  scale_fill_manual(name = "", values = c(MCMC = "red", VBfixphi = "blue",VBfixssqphi="green"))+
  scale_color_manual(name = "",values = c(MCMC = "red", VBfixphi = "blue",VBfixssqphi="green"))+
  geom_vline(xintercept =sigma2, color="black",size=1.5)+
  coord_cartesian(xlim=c(0,2))+
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
xSeqsigma2<-seq(0.1,6, length.out=length(newdensitySigma2))
set.seed(5153527)
foosigma2<-sample(x_sigma2,size = iter.thetaphissqphi, replace = T, prob = ELBOsigma2)

foosigma2<-data.frame(foosigma2)
foosigma2fixphi<-data.frame(xSeqsigma2,newdensitySigma2)
MCMC_sigma2<-parameterMatrix[-(1:burninperiod),4]
MCMC_sigma2<-data.frame(MCMC_sigma2)
colnames(foosigma2)[1]<-'VBfixssqphi'
colnames(foosigma2fixphi)[1]<-'xSeqsigma2';colnames(foosigma2fixphi)[2]<-'newdensitySigma2'
colnames(MCMC_sigma2)[1]<-'MCMC_sigma2'
head(foosigma2)
head(foosigma2fixphi)
head(MCMC_sigma2)

ggplot()+
  geom_density(data=MCMC_sigma2, aes(x=MCMC_sigma2, fill="MCMC",color="MCMC"),alpha=0.2, linetype="solid",size=1.5)+
  geom_line(data=foosigma2fixphi, aes(x=xSeqsigma2,y=newdensitySigma2,fill="VBfixphi",color="VBfixphi"),alpha=0.6,linetype="solid",size=1.5)+
  geom_density(data=foosigma2, aes(x=VBfixssqphi,fill="VBfixssqphi",color="VBfixssqphi"),alpha=0.4,linetype="solid",size=1.5)+
  scale_fill_manual(name = "", values = c(MCMC = "red", VBfixphi = "blue",VBfixssqphi="green"))+
  scale_color_manual(name = "",values = c(MCMC = "red", VBfixphi = "blue",VBfixssqphi="green"))+
  geom_vline(xintercept =sigma2, color="black",size=1)+
  coord_cartesian(xlim=c(0,2))+
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
## Phi
################################################################################
foophi<-sample(thetaphissqphi,size = iter.thetaphissqphi, replace = T, prob = ELBOphi )
set.seed(5153527)
foo<-sample(thetaphi,size = iter.thetaphi, replace = T, prob = ELBOvectorweights )

foo<-data.frame(foo) #VBfixphi
foophi<-data.frame(foophi) # VBfixssqphi
MCMC_phi<-parameterMatrix[-(1:burninperiod),3] #MCMCphi
MCMC_phi<-data.frame(MCMC_phi)

colnames(foo)[1]<-'VBfixphi'
colnames(foophi)[1]<-'VBfixssqphi'
colnames(MCMC_phi)[1]<-'MCMCphi'

ggplot()+
  geom_density(data=MCMC_phi, aes(x=MCMCphi, fill="MCMC",color="MCMC"),alpha=0.2, linetype="solid",size=1.5)+
  geom_density(data=foo, aes(x=VBfixphi,fill="VBfixphi",color="VBfixphi"),alpha=0.3,linetype="solid",size=1.5)+
  geom_density(data=foophi, aes(x=VBfixssqphi, fill="VBfixssqphi",color="VBfixssqphi"),alpha=0.4,linetype="solid",size=1.5)+
  scale_fill_manual(name = "", values = c(MCMC = "red", VBfixphi = "blue",VBfixssqphi="green"))+
  scale_color_manual(name = "",values = c(MCMC = "red", VBfixphi = "blue",VBfixssqphi="green")) +
  geom_vline(xintercept =phi, color="black",size=2)+
  coord_cartesian(xlim=c(0,1.4))+
  theme(legend.position = "right")+
  ggtitle("Phi")+
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
## W3 
################################################################################
xSeqW3<-seq(-3,3, length.out=iter.thetaphi)
xSeqW3ssqphi<-seq(-3,3, length.out=length(newdensityW3ssqphi))


W3mcmc<-wMat[-(1:burninperiod),3]
df_W3_1<-data.frame(W3mcmc)
df_W3_2<-data.frame(xSeqW3,newdensityW3)
df_W3_3<-data.frame(xSeqW3ssqphi,newdensityW3ssqphi)
colnames(df_W3_1)[1]<-'MCMC_W3'


ggplot()+
  geom_density(data=df_W3_1, aes(x=MCMC_W3, fill="MCMC",color="MCMC"),alpha=0.2, linetype="solid",size=1.5)+
  geom_line(data=df_W3_2, aes(x=xSeqW3,y=newdensityW3,fill="VBfixphi",color="VBfixphi"),alpha=0.6,linetype="solid",size=1.5)+
  geom_line(data=df_W3_3, aes(x=xSeqW3ssqphi,y=newdensityW3ssqphi,fill="VBfixssqphi",color="VBfixssqphi"),alpha=0.6,linetype="solid",size=1.5)+
  scale_fill_manual(name = "", values = c(MCMC = "red", VBfixphi = "blue",VBfixssqphi="green"))+
  scale_color_manual(name = "",values = c(MCMC = "red", VBfixphi = "blue",VBfixssqphi="green"))+
  geom_vline(xintercept =w[3], color="black",size=1.5)+
  coord_cartesian(xlim=c(-3,3))+
  theme(legend.position = "right")+
  ggtitle("W3")+
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
## W25
################################################################################
xSeqW25<-seq(-3,3, length.out=iter.thetaphi)
xSeqW25ssqphi<-seq(-3,3, length.out=length(newdensityW25ssqphi))


W25mcmc<-wMat[-(1:burninperiod),25]
df_W25_1<-data.frame(W25mcmc)
df_W25_2<-data.frame(xSeqW25,newdensityW25)
df_W25_3<-data.frame(xSeqW25ssqphi,newdensityW25ssqphi)
colnames(df_W25_1)[1]<-'MCMC_W25'


ggplot()+
  geom_density(data=df_W25_1, aes(x=MCMC_W25, fill="MCMC",color="MCMC"),alpha=0.2, linetype="solid",size=1.5)+
  geom_line(data=df_W25_2, aes(x=xSeqW25,y=newdensityW25,fill="VBfixphi",color="VBfixphi"),alpha=0.6,linetype="solid",size=1.5)+
  geom_line(data=df_W25_3, aes(x=xSeqW25ssqphi,y=newdensityW25ssqphi,fill="VBfixssqphi",color="VBfixssqphi"),alpha=0.6,linetype="solid",size=1.5)+
  scale_fill_manual(name = "", values = c(MCMC = "red", VBfixphi = "blue",VBfixssqphi="green"))+
  scale_color_manual(name = "",values = c(MCMC = "red", VBfixphi = "blue",VBfixssqphi="green"))+
  geom_vline(xintercept =w[25], color="black",size=1.5)+
  coord_cartesian(xlim=c(-3,3))+
  theme(legend.position = "right")+
  ggtitle("W25")+
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
aucVal
aucVal_VB
aucVal_VBssqphi

ptFinal # Walltime - MCMC
VB_ptFinal # Walltime - INFVB
VB_ptFinalssqphi

ptFinal/VB_ptFinalssqphi # Computational Speedup=(MCMC/INFVB2)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB1)
