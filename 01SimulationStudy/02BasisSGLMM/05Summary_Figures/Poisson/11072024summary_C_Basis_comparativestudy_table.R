################################################################################
################################################################################
# 11072024 for Comparative study Table_Count_Basis
# Stan and INLA should be added later
################################################################################
################################################################################
rm(list=ls())
setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/02BasisSGLMM/05Summary_Figures")
# dev.off()
getwd()
################################################################################
################################################################################
################################################################################
################################################################################
# 20 Basis functions 
################################################################################
################################################################################
# When phi=0.1
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/phi1_Gen25k_1SpatialData.RData")
load(file="../01DataGeneration/phi1_Gen25k_1EigenBasis.RData")
load(file="../02MCMC/Poisson/Basis_C_phi1_Gen25k_1_MCMC.RData")
load(file="../03INFVB_ParallelComputing/Poisson/Basis_C_phi1_Gen25k_1_VBssqParallel.RData")
load(file="../04HybridMFVB/Poisson/Basis_C_phi1_Gen25k_1_HMFVB.RData")
load(file="../06INLA/Poisson/Basis_C_phi1_Gen25k_1_INLA.RData")
load(file="../10Stan/Poisson/Basis_C_phi1_Gen25k_1_MCMC_Stan.RData")


obs_RMSPE # RMSPE - MCMC
obs_RMSPE_MFVB # RMSPE - MFVB
obs_RMSPE_VB # RMSPE - INFVB
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB
TotTime_INLA # Walltime - INLA
totTime # Walltime - Stan

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)
ptFinal/TotTime_INLA
ptFinal/totTime

################################################################################
################################################################################
# When phi=0.3
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/phi2_Gen25k_1SpatialData.RData")
load(file="../01DataGeneration/phi2_Gen25k_1EigenBasis.RData")
load(file="../02MCMC/Poisson/Basis_C_phi2_Gen25k_1_MCMC.RData")
load(file="../03INFVB_ParallelComputing/Poisson/Basis_C_phi2_Gen25k_1_VBssqParallel.RData")
load(file="../04HybridMFVB/Poisson/Basis_C_phi2_Gen25k_1_HMFVB.RData")
load(file="../06INLA/Poisson/Basis_C_phi2_Gen25k_1_INLA.RData")
load(file="../10Stan/Poisson/Basis_C_phi2_Gen25k_1_MCMC_Stan.RData")


obs_RMSPE # RMSPE - MCMC
obs_RMSPE_MFVB # RMSPE - MFVB
obs_RMSPE_VB # RMSPE - INFVB
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB
TotTime_INLA # Walltime - INLA
totTime # Walltime - Stan

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)
ptFinal/TotTime_INLA
ptFinal/totTime

################################################################################
################################################################################
# When phi=0.5
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/phi3_Gen25k_1SpatialData.RData")
load(file="../01DataGeneration/phi3_Gen25k_1EigenBasis.RData")
load(file="../02MCMC/Poisson/Basis_C_phi3_Gen25k_1_MCMC.RData")
load(file="../03INFVB_ParallelComputing/Poisson/Basis_C_phi3_Gen25k_1_VBssqParallel.RData")
load(file="../04HybridMFVB/Poisson/Basis_C_phi3_Gen25k_1_HMFVB.RData")
load(file="../06INLA/Poisson/Basis_C_phi3_Gen25k_1_INLA.RData")
load(file="../10Stan/Poisson/Basis_C_phi3_Gen25k_1_MCMC_Stan.RData")


obs_RMSPE # RMSPE - MCMC
obs_RMSPE_MFVB # RMSPE - MFVB
obs_RMSPE_VB # RMSPE - INFVB
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB
TotTime_INLA # Walltime - INLA
totTime # Walltime - Stan

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)
ptFinal/TotTime_INLA
ptFinal/totTime
################################################################################
################################################################################
# When phi=0.7
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/phi4_Gen25k_1SpatialData.RData")
load(file="../01DataGeneration/phi4_Gen25k_1EigenBasis.RData")
load(file="../02MCMC/Poisson/Basis_C_phi4_Gen25k_1_MCMC.RData")
load(file="../03INFVB_ParallelComputing/Poisson/Basis_C_phi4_Gen25k_1_VBssqParallel.RData")
load(file="../04HybridMFVB/Poisson/Basis_C_phi4_Gen25k_1_HMFVB.RData")
load(file="../06INLA/Poisson/Basis_C_phi4_Gen25k_1_INLA.RData")
load(file="../10Stan/Poisson/Basis_C_phi4_Gen25k_1_MCMC_Stan.RData")


obs_RMSPE # RMSPE - MCMC
obs_RMSPE_MFVB # RMSPE - MFVB
obs_RMSPE_VB # RMSPE - INFVB
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB
TotTime_INLA # Walltime - INLA
totTime # Walltime - Stan

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)
ptFinal/TotTime_INLA
ptFinal/totTime





################################################################################
################################################################################
################################################################################
################################################################################
# 50 Basis functions 
################################################################################
################################################################################


################################################################################
# When phi=0.1
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/phi1_Gen25k_1SpatialData.RData")
load(file="../01DataGeneration/phi1_Gen25k_1EigenBasis.RData")
load(file="../02MCMC/Poisson/50Basis_C_phi1_Gen25k_1_MCMC.RData")
load(file="../03INFVB_ParallelComputing/Poisson/50Basis_C_phi1_Gen25k_1_VBssqParallel.RData")
load(file="../04HybridMFVB/Poisson/50Basis_C_phi1_Gen25k_1_HMFVB.RData")
load(file="../06INLA/Poisson/50Basis_C_phi1_Gen25k_1_INLA.RData")
load(file="../10Stan/Poisson/50Basis_C_phi1_Gen25k_1_MCMC_Stan.RData")


obs_RMSPE # RMSPE - MCMC
obs_RMSPE_MFVB # RMSPE - MFVB
obs_RMSPE_VB # RMSPE - INFVB
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB
TotTime_INLA # Walltime - INLA
totTime # Walltime - Stan

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)
ptFinal/TotTime_INLA
ptFinal/totTime

################################################################################
################################################################################
# When phi=0.3
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/phi2_Gen25k_1SpatialData.RData")
load(file="../01DataGeneration/phi2_Gen25k_1EigenBasis.RData")
load(file="../02MCMC/Poisson/50Basis_C_phi2_Gen25k_1_MCMC.RData")
load(file="../03INFVB_ParallelComputing/Poisson/50Basis_C_phi2_Gen25k_1_VBssqParallel.RData")
load(file="../04HybridMFVB/Poisson/50Basis_C_phi2_Gen25k_1_HMFVB.RData")
load(file="../06INLA/Poisson/50Basis_C_phi2_Gen25k_1_INLA.RData")
load(file="../10Stan/Poisson/50Basis_C_phi2_Gen25k_1_MCMC_Stan.RData")


obs_RMSPE # RMSPE - MCMC
obs_RMSPE_MFVB # RMSPE - MFVB
obs_RMSPE_VB # RMSPE - INFVB
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB
TotTime_INLA # Walltime - INLA
totTime # Walltime - Stan

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)
ptFinal/TotTime_INLA
ptFinal/totTime

################################################################################
################################################################################
# When phi=0.5
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/phi3_Gen25k_1SpatialData.RData")
load(file="../01DataGeneration/phi3_Gen25k_1EigenBasis.RData")
load(file="../02MCMC/Poisson/50Basis_C_phi3_Gen25k_1_MCMC.RData")
load(file="../03INFVB_ParallelComputing/Poisson/50Basis_C_phi3_Gen25k_1_VBssqParallel.RData")
load(file="../04HybridMFVB/Poisson/50Basis_C_phi3_Gen25k_1_HMFVB.RData")
load(file="../06INLA/Poisson/50Basis_C_phi3_Gen25k_1_INLA.RData")
load(file="../10Stan/Poisson/50Basis_C_phi3_Gen25k_1_MCMC_Stan.RData")


obs_RMSPE # RMSPE - MCMC
obs_RMSPE_MFVB # RMSPE - MFVB
obs_RMSPE_VB # RMSPE - INFVB
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB
TotTime_INLA # Walltime - INLA
totTime # Walltime - Stan

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)
ptFinal/TotTime_INLA
ptFinal/totTime

################################################################################
################################################################################
# When phi=0.7
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/phi4_Gen25k_1SpatialData.RData")
load(file="../01DataGeneration/phi4_Gen25k_1EigenBasis.RData")
load(file="../02MCMC/Poisson/50Basis_C_phi4_Gen25k_1_MCMC.RData")
load(file="../03INFVB_ParallelComputing/Poisson/50Basis_C_phi4_Gen25k_1_VBssqParallel.RData")
load(file="../04HybridMFVB/Poisson/50Basis_C_phi4_Gen25k_1_HMFVB.RData")
load(file="../06INLA/Poisson/50Basis_C_phi4_Gen25k_1_INLA.RData")
load(file="../10Stan/Poisson/50Basis_C_phi4_Gen25k_1_MCMC_Stan.RData")


obs_RMSPE # RMSPE - MCMC
obs_RMSPE_MFVB # RMSPE - MFVB
obs_RMSPE_VB # RMSPE - INFVB
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB
TotTime_INLA # Walltime - INLA
totTime # Walltime - Stan

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)
ptFinal/TotTime_INLA
ptFinal/totTime










################################################################################
################################################################################
################################################################################
################################################################################
# 100 Basis functions 
################################################################################
################################################################################



################################################################################
# When phi=0.1
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/phi1_Gen25k_1SpatialData.RData")
load(file="../01DataGeneration/phi1_Gen25k_1EigenBasis.RData")
load(file="../02MCMC/Poisson/100Basis_C_phi1_Gen25k_1_MCMC.RData")
load(file="../03INFVB_ParallelComputing/Poisson/100Basis_C_phi1_Gen25k_1_VBssqParallel.RData")
load(file="../04HybridMFVB/Poisson/100Basis_C_phi1_Gen25k_1_HMFVB.RData")
load(file="../06INLA/Poisson/100Basis_C_phi1_Gen25k_1_INLA.RData")
load(file="../10Stan/Poisson/100Basis_C_phi1_Gen25k_1_MCMC_Stan.RData")


obs_RMSPE # RMSPE - MCMC
obs_RMSPE_MFVB # RMSPE - MFVB
obs_RMSPE_VB # RMSPE - INFVB
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB
TotTime_INLA # Walltime - INLA
totTime # Walltime - Stan

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)
ptFinal/TotTime_INLA
ptFinal/totTime


################################################################################
################################################################################
# When phi=0.3
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/phi2_Gen25k_1SpatialData.RData")
load(file="../01DataGeneration/phi2_Gen25k_1EigenBasis.RData")
load(file="../02MCMC/Poisson/100Basis_C_phi2_Gen25k_1_MCMC.RData")
load(file="../03INFVB_ParallelComputing/Poisson/100Basis_C_phi2_Gen25k_1_VBssqParallel.RData")
load(file="../04HybridMFVB/Poisson/100Basis_C_phi2_Gen25k_1_HMFVB.RData")
load(file="../06INLA/Poisson/100Basis_C_phi2_Gen25k_1_INLA.RData")
load(file="../10Stan/Poisson/100Basis_C_phi2_Gen25k_1_MCMC_Stan.RData")


obs_RMSPE # RMSPE - MCMC
obs_RMSPE_MFVB # RMSPE - MFVB
obs_RMSPE_VB # RMSPE - INFVB
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB
TotTime_INLA # Walltime - INLA
totTime # Walltime - Stan

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)
ptFinal/TotTime_INLA
ptFinal/totTime

################################################################################
################################################################################
# When phi=0.5
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/phi3_Gen25k_1SpatialData.RData")
load(file="../01DataGeneration/phi3_Gen25k_1EigenBasis.RData")
load(file="../02MCMC/Poisson/100Basis_C_phi3_Gen25k_1_MCMC.RData")
load(file="../03INFVB_ParallelComputing/Poisson/100Basis_C_phi3_Gen25k_1_VBssqParallel.RData")
load(file="../04HybridMFVB/Poisson/100Basis_C_phi3_Gen25k_1_HMFVB.RData")
load(file="../06INLA/Poisson/100Basis_C_phi3_Gen25k_1_INLA.RData")
load(file="../10Stan/Poisson/100Basis_C_phi3_Gen25k_1_MCMC_Stan.RData")


obs_RMSPE # RMSPE - MCMC
obs_RMSPE_MFVB # RMSPE - MFVB
obs_RMSPE_VB # RMSPE - INFVB
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB
TotTime_INLA # Walltime - INLA
totTime # Walltime - Stan

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)
ptFinal/TotTime_INLA
ptFinal/totTime



################################################################################
################################################################################
# When phi=0.7
################################################################################
################################################################################
rm(list=ls())
load(file="../01DataGeneration/phi4_Gen25k_1SpatialData.RData")
load(file="../01DataGeneration/phi4_Gen25k_1EigenBasis.RData")
load(file="../02MCMC/Poisson/100Basis_C_phi4_Gen25k_1_MCMC.RData")
load(file="../03INFVB_ParallelComputing/Poisson/100Basis_C_phi4_Gen25k_1_VBssqParallel.RData")
load(file="../04HybridMFVB/Poisson/100Basis_C_phi4_Gen25k_1_HMFVB.RData")
load(file="../06INLA/Poisson/100Basis_C_phi4_Gen25k_1_INLA.RData")
load(file="../10Stan/Poisson/100Basis_C_phi4_Gen25k_1_MCMC_Stan.RData")


obs_RMSPE # RMSPE - MCMC
obs_RMSPE_MFVB # RMSPE - MFVB
obs_RMSPE_VB # RMSPE - INFVB
obs_RMSPE_INLA
obs_RMSPE_MCMC_Stan


ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB
TotTime_INLA # Walltime - INLA
totTime # Walltime - Stan

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)
ptFinal/TotTime_INLA
ptFinal/totTime

