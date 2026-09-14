################################################################################
################################################################################
# Summary Gaussian for CImethod
################################################################################
################################################################################
rm(list=ls())
getwd()
# setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/03Normal_Supplementary/07Measures")
################################################################################
################################################################################
#phi1 50 Simulations
load("G_phi1_Gen500_TVCI.RData")

length(TVCItest_MCMC_mean_vector)
length(TVCItest_VBphi_mean_vector)

mean(TVCItest_MCMC_mean_vector)
mean(TVCItest_VBphi_mean_vector)

################################################################################
################################################################################
#phi2 50 Simulations
load("G_phi2_Gen500_TVCI.RData")
length(TVCItest_MCMC_mean_vector)
length(TVCItest_VBphi_mean_vector)

mean(TVCItest_MCMC_mean_vector)
mean(TVCItest_VBphi_mean_vector)


################################################################################
################################################################################
#phi3 50 Simulations
load("G_phi3_Gen500_TVCI.RData")
length(TVCItest_MCMC_mean_vector)
length(TVCItest_VBphi_mean_vector)

mean(TVCItest_MCMC_mean_vector)
mean(TVCItest_VBphi_mean_vector)



################################################################################
################################################################################
#phi4 50 Simulations
load("G_phi4_Gen500_TVCI.RData")
length(TVCItest_MCMC_mean_vector)
length(TVCItest_VBphi_mean_vector)

mean(TVCItest_MCMC_mean_vector)
mean(TVCItest_VBphi_mean_vector)

