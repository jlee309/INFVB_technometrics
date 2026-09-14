################################################################################
# 11142024 for Table2 with 50 simulations. 
# Binary Table_Basis
################################################################################
################################################################################
rm(list=ls())
getwd()
# setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/02BasisSGLMM/05Summary_Figures/Binary")
################################################################################
################################################################################

n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  
  AUC_MCMCvector<-c()
  AUC_INLAvector<-c()
  AUC_VBssqvector<-c()
  AUC_HMFVBvector<-c()
  
  comptime_MCMCvector<-c()
  comptime_INLAvector<-c()
  comptime_VBssqvector<-c()
  comptime_HMFVBvector<-c()
  
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    
    file_name1 <- paste0("../../01DataGeneration/phi", phiselect, "_Gen25k_", sim, "SpatialData.RData")
    file_name2 <- paste0("../../01DataGeneration/phi", phiselect, "_Gen25k_", sim, "EigenBasis.RData")
    file_name3 <- paste0("../../02MCMC/Binary/Basis_B_phi", phiselect, "_Gen25k_", sim, "_MCMC.RData")
    file_name4 <- paste0("../../03INFVB_ParallelComputing/Binary/Basis_B_phi", phiselect, "_Gen25k_", sim, "_VBssqParallel.RData")
    file_name5 <- paste0("../../04HybridMFVB/Binary/Basis_B_phi", phiselect, "_Gen25k_", sim, "_HMFVB.RData")
    file_name6 <- paste0("../../06INLA/Binary/Basis_B_phi", phiselect, "_Gen25k_", sim, "_INLA.RData")    
    load(file_name1)
    load(file_name2)
    load(file_name3)
    load(file_name4)
    load(file_name5)
    load(file_name6)
    
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    
    AUC_MCMC<-aucVal
    AUC_VBssq<-aucVal_VB
    AUC_INLA<-aucVal_INLA
    AUC_HMFVB<-aucVal_MFVB
    
    comptime_MCMC<- ptFinal[3] # Walltime - MCMC
    comptime_VBssq<-VB_ptFinal[3] # Walltime - INFVB_ssq
    comptime_INLA<-TotTime_INLA # Walltime - INLA
    comptime_HMFVB<-MFVB_ptFinal[3]
    
    AUC_MCMCvector<-c(AUC_MCMCvector,AUC_MCMC)
    AUC_INLAvector<-c(AUC_INLAvector,AUC_INLA)
    AUC_VBssqvector<-c(AUC_VBssqvector, AUC_VBssq)
    AUC_HMFVBvector<-c(AUC_HMFVBvector,AUC_HMFVB)
    
    comptime_MCMCvector<-c(comptime_MCMCvector,comptime_MCMC)
    comptime_INLAvector<-c(comptime_INLAvector,comptime_INLA)
    comptime_VBssqvector<-c(comptime_VBssqvector, comptime_VBssq)
    comptime_HMFVBvector<-c(comptime_HMFVBvector,comptime_HMFVB)
  }
  save(AUC_MCMCvector, AUC_INLAvector, 
       AUC_VBssqvector, AUC_HMFVBvector,
       comptime_MCMCvector, comptime_INLAvector, 
       comptime_VBssqvector, comptime_HMFVBvector,
       file = paste0("Basis_B_phi",phiselect,"_Gen25k_nsim50_Table.RData"))
}

