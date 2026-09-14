################################################################################
################################################################################
# Summary CRPS for Binary Basis
################################################################################
################################################################################
rm(list=ls())
getwd()
# setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/02BasisSGLMM/11Measures/Binary")
################################################################################
################################################################################
#50 Basis 
################################################################################
################################################################################
#phi1 50 Simulations
rm(list=ls())
load("50Basis_B_phi1_Gen25k_CRPS.RData")
length(crpsTotal_MCMC_mean_vector)
length(crpsTotal_HMFVB_mean_vector)
length(crpsTotal_VBssq_mean_vector)
length(crpsTotal_INLA_mean_vector)


crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_HMFVB_mean_final<-mean(crpsTotal_HMFVB_mean_vector)
crpsTotal_VBssq_mean_final<-mean(crpsTotal_VBssq_mean_vector)
crpsTotal_INLA_mean_final<-mean(crpsTotal_INLA_mean_vector)

logsTotal_MCMC_median_final<-median(logsTotal_MCMC_median_vector)
logsTotal_VBssq_median_final<-median(logsTotal_VBssq_median_vector)
logsTotal_HMFVB_median_final<-median(logsTotal_HMFVB_median_vector)
logsTotal_INLA_meadian_final<-median(logsTotal_INLA_meadian_vector)

crpsTotal_MCMC_mean_final
crpsTotal_HMFVB_mean_final
crpsTotal_VBssq_mean_final
crpsTotal_INLA_mean_final

# logsTotal_MCMC_median_final
# logsTotal_VBssq_median_final
# logsTotal_HMFVB_median_final
# logsTotal_INLA_meadian_final


################################################################################
################################################################################
#50 Basis 
################################################################################
################################################################################
#phi2 50 Simulations
load("50Basis_B_phi2_Gen25k_CRPS.RData")
length(crpsTotal_MCMC_mean_vector)
length(crpsTotal_HMFVB_mean_vector)
length(crpsTotal_VBssq_mean_vector)
length(crpsTotal_INLA_mean_vector)

crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_HMFVB_mean_final<-mean(crpsTotal_HMFVB_mean_vector)
crpsTotal_VBssq_mean_final<-mean(crpsTotal_VBssq_mean_vector)
crpsTotal_INLA_mean_final<-mean(crpsTotal_INLA_mean_vector)

logsTotal_MCMC_median_final<-median(logsTotal_MCMC_median_vector)
logsTotal_VBssq_median_final<-median(logsTotal_VBssq_median_vector)
logsTotal_HMFVB_median_final<-median(logsTotal_HMFVB_median_vector)
logsTotal_INLA_meadian_final<-median(logsTotal_INLA_meadian_vector)

crpsTotal_MCMC_mean_final
crpsTotal_HMFVB_mean_final
crpsTotal_VBssq_mean_final
crpsTotal_INLA_mean_final

# logsTotal_MCMC_median_final
# logsTotal_VBssq_median_final
# logsTotal_HMFVB_median_final
# logsTotal_INLA_meadian_final



################################################################################
################################################################################
#50 Basis 
################################################################################
################################################################################
#phi3 50 Simulations
load("50Basis_B_phi3_Gen25k_CRPS.RData")
length(crpsTotal_MCMC_mean_vector)
length(crpsTotal_HMFVB_mean_vector)
length(crpsTotal_VBssq_mean_vector)
length(crpsTotal_INLA_mean_vector)

crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_HMFVB_mean_final<-mean(crpsTotal_HMFVB_mean_vector)
crpsTotal_VBssq_mean_final<-mean(crpsTotal_VBssq_mean_vector)
crpsTotal_INLA_mean_final<-mean(crpsTotal_INLA_mean_vector)

logsTotal_MCMC_median_final<-median(logsTotal_MCMC_median_vector)
logsTotal_VBssq_median_final<-median(logsTotal_VBssq_median_vector)
logsTotal_HMFVB_median_final<-median(logsTotal_HMFVB_median_vector)
logsTotal_INLA_meadian_final<-median(logsTotal_INLA_meadian_vector)

crpsTotal_MCMC_mean_final
crpsTotal_HMFVB_mean_final
crpsTotal_VBssq_mean_final
crpsTotal_INLA_mean_final

# logsTotal_MCMC_median_final
# logsTotal_VBssq_median_final
# logsTotal_HMFVB_median_final
# logsTotal_INLA_meadian_final



################################################################################
################################################################################
#50 Basis 
################################################################################
################################################################################
#phi4 50 Simulations
load("50Basis_B_phi4_Gen25k_CRPS.RData")
length(crpsTotal_MCMC_mean_vector)
length(crpsTotal_HMFVB_mean_vector)
length(crpsTotal_VBssq_mean_vector)
length(crpsTotal_INLA_mean_vector)

crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_HMFVB_mean_final<-mean(crpsTotal_HMFVB_mean_vector)
crpsTotal_VBssq_mean_final<-mean(crpsTotal_VBssq_mean_vector)
crpsTotal_INLA_mean_final<-mean(crpsTotal_INLA_mean_vector)

logsTotal_MCMC_median_final<-median(logsTotal_MCMC_median_vector)
logsTotal_VBssq_median_final<-median(logsTotal_VBssq_median_vector)
logsTotal_HMFVB_median_final<-median(logsTotal_HMFVB_median_vector)
logsTotal_INLA_meadian_final<-median(logsTotal_INLA_meadian_vector)

crpsTotal_MCMC_mean_final
crpsTotal_HMFVB_mean_final
crpsTotal_VBssq_mean_final
crpsTotal_INLA_mean_final

# logsTotal_MCMC_median_final
# logsTotal_VBssq_median_final
# logsTotal_HMFVB_median_final
# logsTotal_INLA_meadian_final















################################################################################
################################################################################
#20 Basis 
################################################################################
################################################################################
#phi1 50 Simulations
rm(list=ls())
load("Basis_B_phi1_Gen25k_CRPS.RData")
length(crpsTotal_MCMC_mean_vector)
length(crpsTotal_HMFVB_mean_vector)
length(crpsTotal_VBssq_mean_vector)
length(crpsTotal_INLA_mean_vector)


crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_HMFVB_mean_final<-mean(crpsTotal_HMFVB_mean_vector)
crpsTotal_VBssq_mean_final<-mean(crpsTotal_VBssq_mean_vector)
crpsTotal_INLA_mean_final<-mean(crpsTotal_INLA_mean_vector)

logsTotal_MCMC_median_final<-median(logsTotal_MCMC_median_vector)
logsTotal_VBssq_median_final<-median(logsTotal_VBssq_median_vector)
logsTotal_HMFVB_median_final<-median(logsTotal_HMFVB_median_vector)
logsTotal_INLA_meadian_final<-median(logsTotal_INLA_meadian_vector)

crpsTotal_MCMC_mean_final
crpsTotal_HMFVB_mean_final
crpsTotal_VBssq_mean_final
crpsTotal_INLA_mean_final

# logsTotal_MCMC_median_final
# logsTotal_VBssq_median_final
# logsTotal_HMFVB_median_final
# logsTotal_INLA_meadian_final


################################################################################
################################################################################
#20 Basis 
################################################################################
################################################################################
#phi2 50 Simulations
load("Basis_B_phi2_Gen25k_CRPS.RData")
length(crpsTotal_MCMC_mean_vector)
length(crpsTotal_HMFVB_mean_vector)
length(crpsTotal_VBssq_mean_vector)
length(crpsTotal_INLA_mean_vector)

crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_HMFVB_mean_final<-mean(crpsTotal_HMFVB_mean_vector)
crpsTotal_VBssq_mean_final<-mean(crpsTotal_VBssq_mean_vector)
crpsTotal_INLA_mean_final<-mean(crpsTotal_INLA_mean_vector)

logsTotal_MCMC_median_final<-median(logsTotal_MCMC_median_vector)
logsTotal_VBssq_median_final<-median(logsTotal_VBssq_median_vector)
logsTotal_HMFVB_median_final<-median(logsTotal_HMFVB_median_vector)
logsTotal_INLA_meadian_final<-median(logsTotal_INLA_meadian_vector)

crpsTotal_MCMC_mean_final
crpsTotal_HMFVB_mean_final
crpsTotal_VBssq_mean_final
crpsTotal_INLA_mean_final

# logsTotal_MCMC_median_final
# logsTotal_VBssq_median_final
# logsTotal_HMFVB_median_final
# logsTotal_INLA_meadian_final



################################################################################
################################################################################
#20 Basis 
################################################################################
################################################################################
#phi3 50 Simulations
rm(list=ls())
load("Basis_B_phi3_Gen25k_CRPS.RData")
length(crpsTotal_MCMC_mean_vector)
length(crpsTotal_HMFVB_mean_vector)
length(crpsTotal_VBssq_mean_vector)
length(crpsTotal_INLA_mean_vector)

crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_HMFVB_mean_final<-mean(crpsTotal_HMFVB_mean_vector)
crpsTotal_VBssq_mean_final<-mean(crpsTotal_VBssq_mean_vector)
crpsTotal_INLA_mean_final<-mean(crpsTotal_INLA_mean_vector)

logsTotal_MCMC_median_final<-median(logsTotal_MCMC_median_vector)
logsTotal_VBssq_median_final<-median(logsTotal_VBssq_median_vector)
logsTotal_HMFVB_median_final<-median(logsTotal_HMFVB_median_vector)
logsTotal_INLA_meadian_final<-median(logsTotal_INLA_meadian_vector)

crpsTotal_MCMC_mean_final
crpsTotal_HMFVB_mean_final
crpsTotal_VBssq_mean_final
crpsTotal_INLA_mean_final

# logsTotal_MCMC_median_final
# logsTotal_VBssq_median_final
# logsTotal_HMFVB_median_final
# logsTotal_INLA_meadian_final



################################################################################
################################################################################
# 20 Basis 
################################################################################
################################################################################
#phi4 50 Simulations
load("Basis_B_phi4_Gen25k_CRPS.RData")
length(crpsTotal_MCMC_mean_vector)
length(crpsTotal_HMFVB_mean_vector)
length(crpsTotal_VBssq_mean_vector)
length(crpsTotal_INLA_mean_vector)

crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_HMFVB_mean_final<-mean(crpsTotal_HMFVB_mean_vector)
crpsTotal_VBssq_mean_final<-mean(crpsTotal_VBssq_mean_vector)
crpsTotal_INLA_mean_final<-mean(crpsTotal_INLA_mean_vector)

logsTotal_MCMC_median_final<-median(logsTotal_MCMC_median_vector)
logsTotal_VBssq_median_final<-median(logsTotal_VBssq_median_vector)
logsTotal_HMFVB_median_final<-median(logsTotal_HMFVB_median_vector)
logsTotal_INLA_meadian_final<-median(logsTotal_INLA_meadian_vector)

crpsTotal_MCMC_mean_final
crpsTotal_HMFVB_mean_final
crpsTotal_VBssq_mean_final
crpsTotal_INLA_mean_final

# logsTotal_MCMC_median_final
# logsTotal_VBssq_median_final
# logsTotal_HMFVB_median_final
# logsTotal_INLA_meadian_final












################################################################################
################################################################################
#100 Basis 
################################################################################
################################################################################
#phi1 50 Simulations
rm(list=ls())
load("100Basis_B_phi1_Gen25k_CRPS.RData")
length(crpsTotal_MCMC_mean_vector)
length(crpsTotal_HMFVB_mean_vector)
length(crpsTotal_VBssq_mean_vector)
length(crpsTotal_INLA_mean_vector)


crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_HMFVB_mean_final<-mean(crpsTotal_HMFVB_mean_vector)
crpsTotal_VBssq_mean_final<-mean(crpsTotal_VBssq_mean_vector)
crpsTotal_INLA_mean_final<-mean(crpsTotal_INLA_mean_vector)

logsTotal_MCMC_median_final<-median(logsTotal_MCMC_median_vector)
logsTotal_VBssq_median_final<-median(logsTotal_VBssq_median_vector)
logsTotal_HMFVB_median_final<-median(logsTotal_HMFVB_median_vector)
logsTotal_INLA_meadian_final<-median(logsTotal_INLA_meadian_vector)

crpsTotal_MCMC_mean_final
crpsTotal_HMFVB_mean_final
crpsTotal_VBssq_mean_final
crpsTotal_INLA_mean_final

# logsTotal_MCMC_median_final
# logsTotal_VBssq_median_final
# logsTotal_HMFVB_median_final
# logsTotal_INLA_meadian_final


################################################################################
################################################################################
#100 Basis 
################################################################################
################################################################################
#phi2 50 Simulations
load("100Basis_B_phi2_Gen25k_CRPS.RData")
length(crpsTotal_MCMC_mean_vector)
length(crpsTotal_HMFVB_mean_vector)
length(crpsTotal_VBssq_mean_vector)
length(crpsTotal_INLA_mean_vector)

crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_HMFVB_mean_final<-mean(crpsTotal_HMFVB_mean_vector)
crpsTotal_VBssq_mean_final<-mean(crpsTotal_VBssq_mean_vector)
crpsTotal_INLA_mean_final<-mean(crpsTotal_INLA_mean_vector)

logsTotal_MCMC_median_final<-median(logsTotal_MCMC_median_vector)
logsTotal_VBssq_median_final<-median(logsTotal_VBssq_median_vector)
logsTotal_HMFVB_median_final<-median(logsTotal_HMFVB_median_vector)
logsTotal_INLA_meadian_final<-median(logsTotal_INLA_meadian_vector)

crpsTotal_MCMC_mean_final
crpsTotal_HMFVB_mean_final
crpsTotal_VBssq_mean_final
crpsTotal_INLA_mean_final

# logsTotal_MCMC_median_final
# logsTotal_VBssq_median_final
# logsTotal_HMFVB_median_final
# logsTotal_INLA_meadian_final



################################################################################
################################################################################
#100 Basis 
################################################################################
################################################################################
#phi3 50 Simulations
load("100Basis_B_phi3_Gen25k_CRPS.RData")
length(crpsTotal_MCMC_mean_vector)
length(crpsTotal_HMFVB_mean_vector)
length(crpsTotal_VBssq_mean_vector)
length(crpsTotal_INLA_mean_vector)

crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_HMFVB_mean_final<-mean(crpsTotal_HMFVB_mean_vector)
crpsTotal_VBssq_mean_final<-mean(crpsTotal_VBssq_mean_vector)
crpsTotal_INLA_mean_final<-mean(crpsTotal_INLA_mean_vector)

logsTotal_MCMC_median_final<-median(logsTotal_MCMC_median_vector)
logsTotal_VBssq_median_final<-median(logsTotal_VBssq_median_vector)
logsTotal_HMFVB_median_final<-median(logsTotal_HMFVB_median_vector)
logsTotal_INLA_meadian_final<-median(logsTotal_INLA_meadian_vector)

crpsTotal_MCMC_mean_final
crpsTotal_HMFVB_mean_final
crpsTotal_VBssq_mean_final
crpsTotal_INLA_mean_final

# logsTotal_MCMC_median_final
# logsTotal_VBssq_median_final
# logsTotal_HMFVB_median_final
# logsTotal_INLA_meadian_final



################################################################################
################################################################################
# 100 Basis 
################################################################################
################################################################################
#phi4 50 Simulations
load("100Basis_B_phi4_Gen25k_CRPS.RData")
length(crpsTotal_MCMC_mean_vector)
length(crpsTotal_HMFVB_mean_vector)
length(crpsTotal_VBssq_mean_vector)
length(crpsTotal_INLA_mean_vector)

crpsTotal_MCMC_mean_final<-mean(crpsTotal_MCMC_mean_vector)
crpsTotal_HMFVB_mean_final<-mean(crpsTotal_HMFVB_mean_vector)
crpsTotal_VBssq_mean_final<-mean(crpsTotal_VBssq_mean_vector)
crpsTotal_INLA_mean_final<-mean(crpsTotal_INLA_mean_vector)

logsTotal_MCMC_median_final<-median(logsTotal_MCMC_median_vector)
logsTotal_VBssq_median_final<-median(logsTotal_VBssq_median_vector)
logsTotal_HMFVB_median_final<-median(logsTotal_HMFVB_median_vector)
logsTotal_INLA_meadian_final<-median(logsTotal_INLA_meadian_vector)

crpsTotal_MCMC_mean_final
crpsTotal_HMFVB_mean_final
crpsTotal_VBssq_mean_final
crpsTotal_INLA_mean_final

# logsTotal_MCMC_median_final
# logsTotal_VBssq_median_final
# logsTotal_HMFVB_median_final
# logsTotal_INLA_meadian_final









