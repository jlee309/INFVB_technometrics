################################################################################
################################################################################
# 11072024 for Comparative study Table_Binary
# Stan and INLA should be added later
################################################################################
################################################################################
rm(list=ls())
getwd()

################################################################################
################################################################################
# When phi=0.1
################################################################################
################################################################################
load(file="../01DataGeneration/B_phi1_Gen500_1.RData")
load(file="../02MCMC/B_phi1_Gen500_1_MCMC.RData")
load(file="../03INFVBfixphi_ParallelComputing/B_phi1_Gen500_1_VBphiParallel.RData")
load(file="../04INFVBfixssqphi_ParallelComputing/B_phi1_Gen500_1_VBphissqphiParallel.RData")
load(file="../06INLA/B_phi1_Gen500_1_INLA.RData")
load(file="../10Stan/B_phi1_Gen500_1_MCMC_Stan.RData")
################################################################################
aucVal
aucVal_VB
aucVal_VBssqphi
aucVal_INLA
aucVal_Stan

ptFinal # Walltime - MCMC
VB_ptFinal # Walltime - INFVB_phi
VB_ptFinalssqphi # Walltime - INFVB_ssqphi
comptTime_INLA[4] # Walltime - INLA
totTime # Walltime - Stan

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
load(file="../01DataGeneration/B_phi2_Gen500_1.RData")
load(file="../02MCMC/B_phi2_Gen500_1_MCMC.RData")
load(file="../03INFVBfixphi_ParallelComputing/B_phi2_Gen500_1_VBphiParallel.RData")
load(file="../04INFVBfixssqphi_ParallelComputing/B_phi2_Gen500_1_VBphissqphiParallel.RData") 
load(file="../06INLA/B_phi2_Gen500_1_INLA.RData")
load(file="../10Stan/B_phi2_Gen500_1_MCMC_Stan.RData")
################################################################################
aucVal
aucVal_VB
aucVal_VBssqphi
aucVal_INLA
aucVal_Stan

ptFinal # Walltime - MCMC
VB_ptFinal # Walltime - INFVB_phi
VB_ptFinalssqphi # Walltime - INFVB_ssqphi
comptTime_INLA[4] # Walltime - INLA
totTime # Walltime - Stan

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
load(file="../01DataGeneration/B_phi3_Gen500_1.RData")
load(file="../02MCMC/B_phi3_Gen500_1_MCMC.RData")
load(file="../03INFVBfixphi_ParallelComputing/B_phi3_Gen500_1_VBphiParallel.RData")
load(file="../04INFVBfixssqphi_ParallelComputing/B_phi3_Gen500_1_VBphissqphiParallel.RData") 
load(file="../06INLA/B_phi3_Gen500_1_INLA.RData")
load(file="../10Stan/B_phi3_Gen500_1_MCMC_Stan.RData")
################################################################################
aucVal
aucVal_VB
aucVal_VBssqphi
aucVal_INLA
aucVal_Stan

ptFinal # Walltime - MCMC
VB_ptFinal # Walltime - INFVB_phi
VB_ptFinalssqphi # Walltime - INFVB_ssqphi
comptTime_INLA[4] # Walltime - INLA
totTime # Walltime - Stan

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
load(file="../01DataGeneration/B_phi4_Gen500_1.RData")
load(file="../02MCMC/B_phi4_Gen500_1_MCMC.RData")
load(file="../03INFVBfixphi_ParallelComputing/B_phi4_Gen500_1_VBphiParallel.RData")
load(file="../04INFVBfixssqphi_ParallelComputing/B_phi4_Gen500_1_VBphissqphiParallel.RData") 
load(file="../06INLA/B_phi4_Gen500_1_INLA.RData")
load(file="../10Stan/B_phi4_Gen500_1_MCMC_Stan.RData")
################################################################################
aucVal
aucVal_VB
aucVal_VBssqphi
aucVal_INLA
aucVal_Stan

ptFinal # Walltime - MCMC
VB_ptFinal # Walltime - INFVB_phi
VB_ptFinalssqphi # Walltime - INFVB_ssqphi
comptTime_INLA[4] # Walltime - INLA
totTime # Walltime - Stan

ptFinal/VB_ptFinalssqphi # Computational Speedup=(MCMC/INFVB2)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB1)
ptFinal/comptTime_INLA[4]
ptFinal/totTime