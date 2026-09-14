################################################################################
# 11142024 for Table1 with 50 simulations. 
# Gaussian Table
################################################################################
################################################################################
# setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/03Normal_Supplementary/04Summary_Figures")
rm(list=ls())
getwd()
################################################################################
################################################################################

n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  
  Obs_RMSPE_MCMCvector<-c()
  Obs_RMSPE_VBphivector<-c()
  
  comptime_MCMCvector<-c()
  comptime_VBphivector<-c()
  
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name1 <- paste0("../01DataGeneration/G_phi", phiselect, "_Gen500_", sim, ".RData")
    file_name2 <- paste0("../02MCMC/G_phi", phiselect, "_Gen500_", sim, "_MCMC.RData")
    file_name3 <- paste0("../03INFVBfixphi_ParallelComputing/G_phi", phiselect, "_Gen500_", sim, "_VBphiParallel.RData")
    
    load(file_name1)
    load(file_name2)
    load(file_name3)
    
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    
    Z_RMSPE_MCMC
    Z_RMSPE_VB

    comptime_MCMC<- ptFinal[3] # Walltime - MCMC
    comptime_VBphi<-VB_ptFinal[3] # Walltime - INFVB_phi
    
    # ptFinal/VB_ptFinalssqphi # Computational Speedup=(MCMC/INFVB2)
    # ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB1)
    # ptFinal/comptTime_INLA[4]
    
    Obs_RMSPE_MCMCvector<-c(Obs_RMSPE_MCMCvector,Z_RMSPE_MCMC)
    Obs_RMSPE_VBphivector<-c(Obs_RMSPE_VBphivector, Z_RMSPE_VB)
    
    comptime_MCMCvector<-c(comptime_MCMCvector,comptime_MCMC)
    comptime_VBphivector<-c(comptime_VBphivector, comptime_VBphi)
  }
  save(Obs_RMSPE_MCMCvector, 
       Obs_RMSPE_VBphivector,
       comptime_MCMCvector, 
       comptime_VBphivector, 
       file = paste0("G_phi",phiselect,"_Gen500_Full_nsim50_Table.RData"))
}


