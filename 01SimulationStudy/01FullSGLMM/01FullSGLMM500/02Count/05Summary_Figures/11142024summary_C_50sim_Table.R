################################################################################
# 11142024 for Table1 with 50 simulations. 
# Count Table
################################################################################
################################################################################
# setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/02Count/05Summary_Figures")
rm(list=ls())
getwd()
################################################################################
################################################################################

n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


# for(phiselect in 1:length(phiset) ){
for(phiselect in 2:3 ){
  cat("Running phiselect", phiselect, "\n")
  
  Obs_RMSPE_MCMCvector<-c()
  Obs_RMSPE_INLAvector<-c()
  Obs_RMSPE_VBphivector<-c()
  Obs_RMSPE_VBssqphivector<-c()
  
  Lambda_RMSPE_MCMCvector<-c()
  Lambda_RMSPE_INLAvector<-c()
  Lambda_RMSPE_VBphivector<-c()
  Lambda_RMSPE_VBssqphivector<-c()
  
  comptime_MCMCvector<-c()
  comptime_INLAvector<-c()
  comptime_VBphivector<-c()
  comptime_VBssqphivector<-c()
  
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name1 <- paste0("../01DataGeneration/C_phi", phiselect, "_Gen500_", sim, ".RData")
    file_name2 <- paste0("../02MCMC/C_phi", phiselect, "_Gen500_", sim, "_MCMC.RData")
    file_name3 <- paste0("../06INLA/C_phi", phiselect, "_Gen500_", sim, "_INLA.RData")
    file_name4 <- paste0("../03INFVBfixphi_ParallelComputing/C_phi", phiselect, "_Gen500_", sim, "_VBphiParallel.RData")
    file_name5 <- paste0("../04INFVBfixssqphi_ParallelComputing/C_phi", phiselect, "_Gen500_", sim, "_VBphissqphiParallel.RData")
    
    load(file_name1)
    load(file_name2)
    load(file_name3)
    load(file_name4)
    load(file_name5)
    
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    
    obs_RMSPE_MCMC
    obs_RMSPE_VB
    obs_RMSPE_VBssqphi
    obs_RMSPE_INLA
    
    lambda_RMSPE_MCMC
    lambda_RMSPE_VB
    lambda_RMSPE_VBssqphi
    lambda_RMSPE_INLA
    
    comptime_MCMC<- ptFinal[3] # Walltime - MCMC
    comptime_VBphi<-VB_ptFinal[3] # Walltime - INFVB_phi
    comptime_VBssqphi<-VB_ptFinalssqphi[3] # Walltime - INFVB_ssqphi
    comptime_INLA<-comptTime_INLA[4] # Walltime - INLA
    
    # ptFinal/VB_ptFinalssqphi # Computational Speedup=(MCMC/INFVB2)
    # ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB1)
    # ptFinal/comptTime_INLA[4]
    
    Obs_RMSPE_MCMCvector<-c(Obs_RMSPE_MCMCvector,obs_RMSPE_MCMC)
    Obs_RMSPE_INLAvector<-c(Obs_RMSPE_INLAvector,obs_RMSPE_INLA)
    Obs_RMSPE_VBphivector<-c(Obs_RMSPE_VBphivector, obs_RMSPE_VB)
    Obs_RMSPE_VBssqphivector<-c(Obs_RMSPE_VBssqphivector, obs_RMSPE_VBssqphi)
    
    Lambda_RMSPE_MCMCvector<-c(Lambda_RMSPE_MCMCvector,lambda_RMSPE_MCMC)
    Lambda_RMSPE_INLAvector<-c(Lambda_RMSPE_INLAvector,lambda_RMSPE_INLA)
    Lambda_RMSPE_VBphivector<-c(Lambda_RMSPE_VBphivector,lambda_RMSPE_VB)
    Lambda_RMSPE_VBssqphivector<-c(Lambda_RMSPE_VBssqphivector,lambda_RMSPE_VBssqphi)
    
    comptime_MCMCvector<-c(comptime_MCMCvector,comptime_MCMC)
    comptime_INLAvector<-c(comptime_INLAvector,comptime_INLA)
    comptime_VBphivector<-c(comptime_VBphivector, comptime_VBphi)
    comptime_VBssqphivector<-c(comptime_VBssqphivector,comptime_VBssqphi)
  }
  save(Obs_RMSPE_MCMCvector, Obs_RMSPE_INLAvector, 
       Obs_RMSPE_VBphivector, Obs_RMSPE_VBssqphivector, 
       comptime_MCMCvector, comptime_INLAvector, 
       comptime_VBphivector, comptime_VBssqphivector, 
       Lambda_RMSPE_MCMCvector,Lambda_RMSPE_INLAvector,
       Lambda_RMSPE_VBphivector,Lambda_RMSPE_VBssqphivector,
       file = paste0("C_phi",phiselect,"_Gen500_Full_nsim50_Table.RData"))
}
# getwd()
# 
# a1<-Obs_RMSPE_MCMCvector 
# a2<-Obs_RMSPE_INLAvector 
# a3<-Obs_RMSPE_VBphivector 
# a4<-Obs_RMSPE_VBssqphivector 
# a5<-comptime_MCMCvector 
# a6<-comptime_INLAvector 
# a7<-comptime_VBphivector 
# a8<-comptime_VBssqphivector 
# a9<-Lambda_RMSPE_MCMCvector
# a10<-Lambda_RMSPE_INLAvector
# a11<-Lambda_RMSPE_VBphivector
# a12<-Lambda_RMSPE_VBssqphivector
# 
# 
# 
# Obs_RMSPE_MCMCvector<-c(Obs_RMSPE_MCMCvector,a1)
# Obs_RMSPE_INLAvector<-c(Obs_RMSPE_INLAvector,a2)
# Obs_RMSPE_VBphivector<-c(Obs_RMSPE_VBphivector,a3)
# Obs_RMSPE_VBssqphivector<-c(Obs_RMSPE_VBssqphivector,a4)
# comptime_MCMCvector<-c(comptime_MCMCvector,a5)
# comptime_INLAvector<-c(comptime_INLAvector,a6)
# comptime_VBphivector<-c(comptime_VBphivector,a7)
# comptime_VBssqphivector<-c(comptime_VBssqphivector,a8)
# Lambda_RMSPE_MCMCvector<-c(Lambda_RMSPE_MCMCvector,a9)
# Lambda_RMSPE_INLAvector<-c(Lambda_RMSPE_INLAvector,a10)
# Lambda_RMSPE_VBphivector<-c(Lambda_RMSPE_VBphivector,a11)
# Lambda_RMSPE_VBssqphivector<-c(Lambda_RMSPE_VBssqphivector,a12)
# 
# save(Obs_RMSPE_MCMCvector, Obs_RMSPE_INLAvector, 
#      Obs_RMSPE_VBphivector, Obs_RMSPE_VBssqphivector, 
#      comptime_MCMCvector, comptime_INLAvector, 
#      comptime_VBphivector, comptime_VBssqphivector, 
#      Lambda_RMSPE_MCMCvector,Lambda_RMSPE_INLAvector,
#      Lambda_RMSPE_VBphivector,Lambda_RMSPE_VBssqphivector,
#      file = "zC_phi4_Gen500_Full_nsim50_Table.RData")
