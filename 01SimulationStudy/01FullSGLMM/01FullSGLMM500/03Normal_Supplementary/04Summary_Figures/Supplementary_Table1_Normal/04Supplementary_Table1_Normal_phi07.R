################################################################################
################################################################################
# SUMMARY RESULTS FOR NORMAL when n=500 Supplementary Table 1 phi=0.7
################################################################################
################################################################################
rm(list=ls())

load(file="../../01DataGeneration/N_05phi07_Gen500.Rdata")
load(file="../../02MCMC/N_phi07_Gen500_Results_MCMC.RData")
load(file="../../03INFVBfixphi_ParallelComputing/N_phi07_Gen500_Results_VBphiParallel.RData")

################################################################################
################################################################################
## RMSPE FINAL RESULT Comparison between MCMC vs VB 
################################################################################
################################################################################
Z_RMSPE_MCMC # RMSPE - MCMC
Z_RMSPE_VB # RMSPE - VB fixphi
ptFinal # Walltime - MCMC
VB_ptFinal # Walltime - VB fixphi
ptFinal/VB_ptFinal # Computational Speedup 
