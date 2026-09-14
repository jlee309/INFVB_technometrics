################################################################################
# 11142024 for Table2 with 50 simulations. 
# Count Table_50Basis
################################################################################
################################################################################
rm(list=ls())
getwd()
# setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/02BasisSGLMM/05Summary_Figures/Poisson")
################################################################################
################################################################################

n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  
  OBS_MCMCvector<-c()
  OBS_INLAvector<-c()
  OBS_VBssqvector<-c()
  OBS_HMFVBvector<-c()
  
  comptime_MCMCvector<-c()
  comptime_INLAvector<-c()
  comptime_VBssqvector<-c()
  comptime_HMFVBvector<-c()
  
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    
    file_name1 <- paste0("../../01DataGeneration/phi", phiselect, "_Gen25k_", sim, "SpatialData.RData")
    file_name2 <- paste0("../../01DataGeneration/phi", phiselect, "_Gen25k_", sim, "EigenBasis.RData")
    file_name3 <- paste0("../../02MCMC/Poisson/50Basis_C_phi", phiselect, "_Gen25k_", sim, "_MCMC.RData")
    file_name4 <- paste0("../../03INFVB_ParallelComputing/Poisson/50Basis_C_phi", phiselect, "_Gen25k_", sim, "_VBssqParallel.RData")
    file_name5 <- paste0("../../04HybridMFVB/Poisson/50Basis_C_phi", phiselect, "_Gen25k_", sim, "_HMFVB.RData")
    file_name6 <- paste0("../../06INLA/Poisson/50Basis_C_phi", phiselect, "_Gen25k_", sim, "_INLA.RData")    
    load(file_name1)
    load(file_name2)
    load(file_name3)
    load(file_name4)
    load(file_name5)
    load(file_name6)
    
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    
    OBS_MCMC<-obs_RMSPE
    OBS_VBssq<-obs_RMSPE_VB
    OBS_INLA<-obs_RMSPE_INLA
    OBS_HMFVB<-obs_RMSPE_MFVB
    
    comptime_MCMC<- ptFinal[3] # Walltime - MCMC
    comptime_VBssq<-VB_ptFinal[3] # Walltime - INFVB_ssq
    comptime_INLA<-TotTime_INLA # Walltime - INLA
    comptime_HMFVB<-MFVB_ptFinal[3]
    
    OBS_MCMCvector<-c(OBS_MCMCvector,OBS_MCMC)
    OBS_INLAvector<-c(OBS_INLAvector,OBS_INLA)
    OBS_VBssqvector<-c(OBS_VBssqvector, OBS_VBssq)
    OBS_HMFVBvector<-c(OBS_HMFVBvector,OBS_HMFVB)
    
    comptime_MCMCvector<-c(comptime_MCMCvector,comptime_MCMC)
    comptime_INLAvector<-c(comptime_INLAvector,comptime_INLA)
    comptime_VBssqvector<-c(comptime_VBssqvector, comptime_VBssq)
    comptime_HMFVBvector<-c(comptime_HMFVBvector,comptime_HMFVB)
  }
  save(OBS_MCMCvector, OBS_INLAvector, 
       OBS_VBssqvector, OBS_HMFVBvector,
       comptime_MCMCvector, comptime_INLAvector, 
       comptime_VBssqvector, comptime_HMFVBvector,
       file = paste0("50Basis_C_phi",phiselect,"_Gen25k_nsim50_Table.RData"))
}

