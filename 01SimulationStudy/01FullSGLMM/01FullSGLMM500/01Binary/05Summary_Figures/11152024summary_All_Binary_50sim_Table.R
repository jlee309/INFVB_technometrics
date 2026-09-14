################################################################################
################################################################################
# 11152024summary_All_Binary_50sim_Table
################################################################################
################################################################################
rm(list=ls())
getwd()
# setwd("/scratch/jlee309/01VBpaper/01SimulationStudy/01FullSGLMM/01FullSGLMM500/01Binary/05Summary_Figures")
################################################################################
################################################################################
#phi1 50 Simulations
load("B_phi1_Gen500_Full_nsim50_Table.RData")

length(AUC_MCMCvector)
length(AUC_INLAvector)
length(AUC_VBphivector)
length(AUC_VBssqphivector)
length(comptime_MCMCvector)
length(comptime_INLAvector)
length(comptime_VBphivector)
length(comptime_VBssqphivector)

mean(AUC_MCMCvector)
mean(AUC_VBphivector)
mean(AUC_VBssqphivector)
mean(AUC_INLAvector)

mean(comptime_MCMCvector)
mean(comptime_VBphivector)
mean(comptime_VBssqphivector)
mean(comptime_INLAvector)

################################################################################
################################################################################
#phi2 50 Simulations
load("B_phi2_Gen500_Full_nsim50_Table.RData")

length(AUC_MCMCvector)
length(AUC_INLAvector)
length(AUC_VBphivector)
length(AUC_VBssqphivector)
length(comptime_MCMCvector)
length(comptime_INLAvector)
length(comptime_VBphivector)
length(comptime_VBssqphivector)

mean(AUC_MCMCvector)
mean(AUC_VBphivector)
mean(AUC_VBssqphivector)
mean(AUC_INLAvector)

mean(comptime_MCMCvector)
mean(comptime_VBphivector)
mean(comptime_VBssqphivector)
mean(comptime_INLAvector)


################################################################################
################################################################################
#phi3 50 Simulations
load("B_phi3_Gen500_Full_nsim50_Table.RData")

length(AUC_MCMCvector)
length(AUC_INLAvector)
length(AUC_VBphivector)
length(AUC_VBssqphivector)
length(comptime_MCMCvector)
length(comptime_INLAvector)
length(comptime_VBphivector)
length(comptime_VBssqphivector)

mean(AUC_MCMCvector)
mean(AUC_VBphivector)
mean(AUC_VBssqphivector)
mean(AUC_INLAvector)

mean(comptime_MCMCvector)
mean(comptime_VBphivector)
mean(comptime_VBssqphivector)
mean(comptime_INLAvector)



################################################################################
################################################################################
#phi4 50 Simulations
load("B_phi4_Gen500_Full_nsim50_Table.RData")

length(AUC_MCMCvector)
length(AUC_INLAvector)
length(AUC_VBphivector)
length(AUC_VBssqphivector)
length(comptime_MCMCvector)
length(comptime_INLAvector)
length(comptime_VBphivector)
length(comptime_VBssqphivector)

mean(AUC_MCMCvector)
mean(AUC_VBphivector)
mean(AUC_VBssqphivector)
mean(AUC_INLAvector)

mean(comptime_MCMCvector)
mean(comptime_VBphivector)
mean(comptime_VBssqphivector)
mean(comptime_INLAvector)