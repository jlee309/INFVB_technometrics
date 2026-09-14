################################################################################
################################################################################
# 11072024 for Comparative study Table_Count
# Stan and INLA should be added later
################################################################################
################################################################################
setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/02Count/05Summary_Figures")
rm(list=ls())
getwd()

################################################################################
################################################################################
# When phi=0.1
################################################################################
################################################################################
load(file="../01DataGeneration/C_phi1_Gen500_1.RData")
load(file="../02MCMC/C_phi1_Gen500_1_MCMC.RData")
load(file="../03INFVBfixphi_ParallelComputing/C_phi1_Gen500_1_VBphiParallel.RData")
load(file="../04INFVBfixssqphi_ParallelComputing/C_phi1_Gen500_1_VBphissqphiParallel.RData") 
load(file="../06INLA/C_phi1_Gen500_1_INLA.RData")
load(file="../10Stan/C_phi1_Gen500_1_MCMC_Stan.RData")
################################################################################
lambda_RMSPE_MCMC
lambda_RMSPE_VB
lambda_RMSPE_VBssqphi
lambda_RMSPE_INLA
lambda_RMSPE_Stan

# We use obs_RMSPE!!
obs_RMSPE_MCMC
obs_RMSPE_VB
obs_RMSPE_VBssqphi
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
VB_ptFinal # Walltime - INFVB
VB_ptFinalssqphi
comptTime_INLA[4]
totTime

ptFinal/VB_ptFinalssqphi # Computational Speedup=(MCMC/INFVB2)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB1)
ptFinal/comptTime_INLA[4]
ptFinal/totTime
################################################################################
################################################################################
# When phi=0.3
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/C_phi2_Gen500_1.RData")
load(file="../02MCMC/C_phi2_Gen500_1_MCMC.RData")
load(file="../03INFVBfixphi_ParallelComputing/C_phi2_Gen500_1_VBphiParallel.RData")
load(file="../04INFVBfixssqphi_ParallelComputing/C_phi2_Gen500_1_VBphissqphiParallel.RData") 
load(file="../06INLA/C_phi2_Gen500_1_INLA.RData")
load(file="../10Stan/C_phi2_Gen500_1_MCMC_Stan.RData")
################################################################################
lambda_RMSPE_MCMC
lambda_RMSPE_VB
lambda_RMSPE_VBssqphi
lambda_RMSPE_INLA
lambda_RMSPE_Stan

# We use obs_RMSPE!!
obs_RMSPE_MCMC
obs_RMSPE_VB
obs_RMSPE_VBssqphi
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
VB_ptFinal # Walltime - INFVB
VB_ptFinalssqphi
comptTime_INLA[4]
totTime

ptFinal/VB_ptFinalssqphi # Computational Speedup=(MCMC/INFVB2)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB1)
ptFinal/comptTime_INLA[4]
ptFinal/totTime
################################################################################
################################################################################
# When phi=0.5
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/C_phi3_Gen500_1.RData")
load(file="../02MCMC/C_phi3_Gen500_1_MCMC.RData")
load(file="../03INFVBfixphi_ParallelComputing/C_phi3_Gen500_1_VBphiParallel.RData")
load(file="../04INFVBfixssqphi_ParallelComputing/C_phi3_Gen500_1_VBphissqphiParallel.RData") 
load(file="../06INLA/C_phi3_Gen500_1_INLA.RData")
load(file="../10Stan/C_phi3_Gen500_1_MCMC_Stan.RData")
################################################################################
lambda_RMSPE_MCMC
lambda_RMSPE_VB
lambda_RMSPE_VBssqphi
lambda_RMSPE_INLA
lambda_RMSPE_Stan

# We use obs_RMSPE!!
obs_RMSPE_MCMC
obs_RMSPE_VB
obs_RMSPE_VBssqphi
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
VB_ptFinal # Walltime - INFVB
VB_ptFinalssqphi
comptTime_INLA[4]
totTime

ptFinal/VB_ptFinalssqphi # Computational Speedup=(MCMC/INFVB2)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB1)
ptFinal/comptTime_INLA[4]
ptFinal/totTime
################################################################################
################################################################################
# When phi=0.7
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/C_phi4_Gen500_1.RData")
load(file="../02MCMC/C_phi4_Gen500_1_MCMC.RData")
load(file="../03INFVBfixphi_ParallelComputing/C_phi4_Gen500_1_VBphiParallel.RData")
load(file="../04INFVBfixssqphi_ParallelComputing/C_phi4_Gen500_1_VBphissqphiParallel.RData") 
load(file="../06INLA/C_phi4_Gen500_1_INLA.RData")
load(file="../10Stan/C_phi4_Gen500_1_MCMC_Stan.RData")
################################################################################
lambda_RMSPE_MCMC
lambda_RMSPE_VB
lambda_RMSPE_VBssqphi
lambda_RMSPE_INLA
lambda_RMSPE_Stan

# We use obs_RMSPE!!
obs_RMSPE_MCMC
obs_RMSPE_VB
obs_RMSPE_VBssqphi
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
VB_ptFinal # Walltime - INFVB
VB_ptFinalssqphi
comptTime_INLA[4]
totTime

ptFinal/VB_ptFinalssqphi # Computational Speedup=(MCMC/INFVB2)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB1)
ptFinal/comptTime_INLA[4]
ptFinal/totTime