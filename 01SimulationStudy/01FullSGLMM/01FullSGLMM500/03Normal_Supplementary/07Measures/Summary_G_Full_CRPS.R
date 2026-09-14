################################################################################
################################################################################
# Summary CRPS for Gaussian
################################################################################
################################################################################
rm(list=ls())
getwd()
# setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/03Normal_Supplementary/07Measures")
################################################################################
################################################################################
#phi1 50 Simulations
load("G_phi1_Gen500_CRPS.RData")
crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_VBphi_mean_final<-mean(crpsTotal_VBphi_mean_vector)

crpsTotal_MCMC_mean_final
crpsTotal_VBphi_mean_final


################################################################################
################################################################################
#phi2 50 Simulations
load("G_phi2_Gen500_CRPS.RData")
crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_VBphi_mean_final<-mean(crpsTotal_VBphi_mean_vector)

crpsTotal_MCMC_mean_final
crpsTotal_VBphi_mean_final

################################################################################
################################################################################
#phi3 50 Simulations
load("G_phi3_Gen500_CRPS.RData")
crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_VBphi_mean_final<-mean(crpsTotal_VBphi_mean_vector)

crpsTotal_MCMC_mean_final
crpsTotal_VBphi_mean_final

################################################################################
################################################################################
#phi4 50 Simulations
load("G_phi4_Gen500_CRPS.RData")
crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_VBphi_mean_final<-mean(crpsTotal_VBphi_mean_vector)

crpsTotal_MCMC_mean_final
crpsTotal_VBphi_mean_final
