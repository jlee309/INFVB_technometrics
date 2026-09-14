################################################################################
# 11142024 for Table1 with 50 simulations. 
# Binary Table
################################################################################
################################################################################
rm(list=ls())
getwd()
# setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/01Binary/05Summary_Figures")
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
  AUC_VBphivector<-c()
  AUC_VBssqphivector<-c()
  
  comptime_MCMCvector<-c()
  comptime_INLAvector<-c()
  comptime_VBphivector<-c()
  comptime_VBssqphivector<-c()
  
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name1 <- paste0("../01DataGeneration/B_phi", phiselect, "_Gen500_", sim, ".RData")
    file_name2 <- paste0("../02MCMC/B_phi", phiselect, "_Gen500_", sim, "_MCMC.RData")
    file_name3 <- paste0("../06INLA/B_phi", phiselect, "_Gen500_", sim, "_INLA.RData")
    file_name4 <- paste0("../03INFVBfixphi_ParallelComputing/B_phi", phiselect, "_Gen500_", sim, "_VBphiParallel.RData")
    file_name5 <- paste0("../04INFVBfixssqphi_ParallelComputing/B_phi", phiselect, "_Gen500_", sim, "_VBphissqphiParallel.RData")
    
    load(file_name1)
    load(file_name2)
    load(file_name3)
    load(file_name4)
    load(file_name5)
    
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    
    AUC_MCMC<-aucVal
    AUC_VBphi<-aucVal_VB
    AUC_VBssqphi<-aucVal_VBssqphi
    AUC_INLA<-aucVal_INLA
    
    comptime_MCMC<- ptFinal[3] # Walltime - MCMC
    comptime_VBphi<-VB_ptFinal[3] # Walltime - INFVB_phi
    comptime_VBssqphi<-VB_ptFinalssqphi[3] # Walltime - INFVB_ssqphi
    comptime_INLA<-comptTime_INLA[4] # Walltime - INLA
    
    # ptFinal/VB_ptFinalssqphi # Computational Speedup=(MCMC/INFVB2)
    # ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB1)
    # ptFinal/comptTime_INLA[4]
    
    AUC_MCMCvector<-c(AUC_MCMCvector,AUC_MCMC)
    AUC_INLAvector<-c(AUC_INLAvector,AUC_INLA)
    AUC_VBphivector<-c(AUC_VBphivector, AUC_VBphi)
    AUC_VBssqphivector<-c(AUC_VBssqphivector, AUC_VBssqphi)
    
    comptime_MCMCvector<-c(comptime_MCMCvector,comptime_MCMC)
    comptime_INLAvector<-c(comptime_INLAvector,comptime_INLA)
    comptime_VBphivector<-c(comptime_VBphivector, comptime_VBphi)
    comptime_VBssqphivector<-c(comptime_VBssqphivector,comptime_VBssqphi)
  }
  save(AUC_MCMCvector, AUC_INLAvector, 
       AUC_VBphivector, AUC_VBssqphivector, 
       comptime_MCMCvector, comptime_INLAvector, 
       comptime_VBphivector, comptime_VBssqphivector, 
       file = paste0("B_phi",phiselect,"_Gen500_Full_nsim50_Table.RData"))
}

