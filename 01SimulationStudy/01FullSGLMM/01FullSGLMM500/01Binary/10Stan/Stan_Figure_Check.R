


rm(list=ls())


load(file="../01DataGeneration/B_phi1_Gen500_1.RData")
load(file="../02MCMC/B_phi1_Gen500_1_MCMC.RData")
load(file="../10Stan/B_phi1_Gen500_1_MCMC_Stan.RData")


aucVal
aucVal_Stan


plot(density(parameterMatrix[,1]))
lines(density(parMat[,1]))

plot(density(parameterMatrix[,2]))
lines(density(parMat[,2]))

plot(density(parameterMatrix[,3]))
lines(density(parMat[,3]))

plot(density(parameterMatrix[,4]))
lines(density(parMat[,4]))

plot(density(wMat[,4]))
lines(density(omegasMat[,4]))

plot(density(wMat[,50]))
lines(density(omegasMat[,50]))


# load(file="../03INFVB_ParallelComputing/Binary/Basis_B_phi1_Gen25k_1_VBssqParallel.RData")
# load(file="../04HybridMFVB/Binary/Basis_B_phi1_Gen25k_1_HMFVB.RData")

obs_RMSPE
lambda_RMSPE 
obs_RMSPE_MCMC_Stan
lambda_RMSPE_Stan

ptFinal # Walltime - MCMC
totTime # Walltime - Stan

ptFinal # Walltime - MCMC
MFVB_ptFinal # Walltime - MFVB
VB_ptFinal # Walltime - INFVB
totTime # Walltime - Stan

ptFinal/MFVB_ptFinal # Computational Speedup=(MCMC/MFVB)
ptFinal/VB_ptFinal # Computational Speedup=(MCMC/INFVB)
ptFinal/totTime # Computational Speedup=(MCMC/Stan)


################################################################################
ess = function(outp,imselags=TRUE)
{
  if (imselags) # truncate number of lags based on imse approach
  {
    chainACov <- acf(outp,type="covariance",plot = FALSE)$acf ## USE AUTOCOVARIANCES
    ACovlen <- length(chainACov)
    gammaACov <- chainACov[1:(ACovlen-1)]+chainACov[2:ACovlen]
    
    m <- 1
    currgamma <- gammaACov[1]
    k <- 1
    while ((k<length(gammaACov)) && (gammaACov[k+1]>0) && (gammaACov[k]>=gammaACov[k+1]))
      k <- k +1
    cat("truncated after ",k," lags\n")
    if (k==length(gammaACov)) # added up until the very last computed autocovariance
      cat("WARNING: may need to compute more autocovariances/autocorrelations for ess\n")
    
    chainACorr = acf(outp,type="correlation",plot = FALSE)$acf ## USE AUTOCORRELATIONS
    if (k==1)
      ACtime = 1
    else
      ACtime <- 1 + 2*sum(chainACorr[2:k])  # add autocorrelations up to lag determined by imse
  }
  else
  {
    chainACorr = acf(outp,type="correlation",plot = FALSE)$acf ## USE AUTOCORRELATIONS
    ACtime <- 1 + 2*sum(chainACorr[-c(1)])
  }
  
  return(length(outp)/ACtime)
}

################################################################################


##########################
# Checking ESS
##########################
#MCMC_Metropolis Hastings 
# Beta1, Beta2 and Sigma2
for( i in 1:dim(parameterMatrix)[2]) {
  print(ess(parameterMatrix[,i]))
}

#MCMC_Metropolis Hastings 
# From Delta1 to Delta20 
for( i in 1:dim(deltaMatrix)[2]) {
  print(ess(deltaMatrix[,i]))
}

##########################
#MCMC_HMC_Stan 
# Beta1, Beta2 and Sigma2
for( i in 1:dim(parMat)[2]) {
  print(ess(parMat[,i]))
}

#MCMC_Stan 
# From Delta1 to Delta20 
for( i in 1:dim(deltaMat)[2]) {
  print(ess(deltaMat[,i]))
}

