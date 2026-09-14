################################################################################
################################################################################
# 10142024
# Running Multiple Simulations for Binary Data
# MCMC
# Test args
################################################################################
################################################################################

rm(list=ls())
library(cli);library(viridis);library(fields);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm)
library(invgamma);library(pROC);library(Matrix);library(emulator)
library(argparser);library(tidyverse);

#Distributing our work to different CPU
args <- arg_parser('jin_test') %>%
add_argument('--task', type = 'integer', help = 'Task number.') %>%
parse_args()

phi_test<-c(1,2,3,4) #it corresponds to c(0.1, 0.3, 0.5, 0.7)
nsim_test<-1:50 #should change this to 50 later
settings<-expand.grid(phi_test=phi_test,nsim_test=nsim_test)
phiselect<-settings$phi_test[args$task]
sim<-settings$nsim_test[args$task]


file_name <- paste0("../01DataGeneration/B_phi", phiselect, "_Gen500_", sim, ".RData")
load(file_name)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed
randomseed <- base_seed + sim + phiselect * 1000
set.seed(randomseed)

    ## Gibbs Sampling 
    n.iteration <- 100000
    mu_beta <- c(0, 0); sigma_beta <- diag(2); alpha1 <- 0.1; beta1 <- 0.1;
    parameterMatrix <- matrix(NA, nrow = n.iteration, ncol = 4)
    parameterMatrix[1,] <- c(1, 1, 1, 0.5)
    inv_sigma_beta <- solve(sigma_beta)
    
    propSD <- 0.1
    lb <- 0
    ub <- sqrt(2)
    
    wMat <- matrix(NA, nrow = n.iteration, ncol = (0.8 * n))
    wMat[1,] <- c(rep(0, (0.8 * n)))
    
    # Training data in R_phi 
    cholR_trainingdata <- chol(R_phi_trainingdata) # Cholesky Decomposition
    invR_trainingdata <- chol2inv(cholR_trainingdata)
    detR_trainingdata <- (-sum(log(diag(cholR_trainingdata)))) # equals -1/2*log(det(R_phi_trainingdata))
    
    pt <- proc.time()
    for (i in 2:n.iteration) {
      if (i %% 100 == 0) { print(i) }
      
      # Estimate W (Every single component)
      wMat[i,] <- wMat[i-1,]
      for (k in 1:(0.8 * n)) {
        proposed_w <- rnorm(1, mean = wMat[i-1,k], sd = 0.2)
        wMat[i,k] <- proposed_w
        num <- (sum(dbinom(x = obs_trainingdata, size = 1, prob = exp(X_trainingdata %*% parameterMatrix[i-1, c(1,2)] + wMat[i,]) / (1 + exp(X_trainingdata %*% parameterMatrix[i-1, c(1,2)] + wMat[i,])), log = TRUE))) - (1/(2*parameterMatrix[i-1,3])) * t(wMat[i,]) %*% (invR_trainingdata %*% wMat[i,])
        wMat[i,k] <- wMat[i-1,k]
        den <- (sum(dbinom(x = obs_trainingdata, size = 1, prob = exp(X_trainingdata %*% parameterMatrix[i-1, c(1,2)] + wMat[i,]) / (1 + exp(X_trainingdata %*% parameterMatrix[i-1, c(1,2)] + wMat[i,])), log = TRUE))) - (1/(2*parameterMatrix[i-1,3])) * t(wMat[i,]) %*% (invR_trainingdata %*% wMat[i,])
        alpha <- min(0, num - den)
        if (alpha > log(runif(1))) {
          wMat[i,k] <- proposed_w
        } else {
          wMat[i,k] <- wMat[i-1,k]
        }
      }
      
      # Estimating Beta
      proposed_beta <- mvrnorm(n = 1, mu = parameterMatrix[i-1, c(1,2)], Sigma = 0.2 * diag(2))
      num <- (sum(dbinom(x = obs_trainingdata, size = 1, prob = exp(X_trainingdata %*% proposed_beta + wMat[i,]) / (1 + exp(X_trainingdata %*% proposed_beta + wMat[i,])), log = TRUE))) - 1/2 * t(proposed_beta) %*% inv_sigma_beta %*% proposed_beta
      den <- (sum(dbinom(x = obs_trainingdata, size = 1, prob = exp(X_trainingdata %*% parameterMatrix[i-1, c(1,2)] + wMat[i,]) / (1 + exp(X_trainingdata %*% parameterMatrix[i-1, c(1,2)] + wMat[i,])), log = TRUE))) - 1/2 * t(parameterMatrix[i-1, c(1,2)]) %*% inv_sigma_beta %*% parameterMatrix[i-1, c(1,2)]
      alpha <- min(0, num - den)
      uniff1 <- runif(1)
      if (alpha > log(uniff1)) {
        parameterMatrix[i, c(1,2)] <- proposed_beta
      } else {
        parameterMatrix[i, c(1,2)] <- parameterMatrix[i-1, c(1,2)]
      }
      
      # Update Sigma^2
      alpha1_tilde <- alpha1 + (0.8 * n) / 2
      beta1_tilde <- beta1 + 1/2 * as.numeric(t(wMat[i,]) %*% invR_trainingdata %*% (wMat[i,]))
      parameterMatrix[i, 3] <- 1 / rgamma(n = 1, shape = alpha1_tilde, rate = beta1_tilde)
      
      # Estimate Phi
      proposed_phi <- rnorm(1, mean = parameterMatrix[i-1, 4], sd = propSD)
      if (proposed_phi < lb | proposed_phi > ub) {
        alpha <- (-1e10)
      } else {
        R_prop_phi <- exp(-distMat[c(1:(0.8 * n)), c(1:(0.8 * n))] / proposed_phi)
        chol_R_prop_phi <- chol(R_prop_phi)
        invR_prop_phi <- chol2inv(chol_R_prop_phi)
        det_R_prop_phi <- (-sum(log(diag(chol_R_prop_phi))))
        num <- det_R_prop_phi - (1/(2 * parameterMatrix[i, 3])) * t(wMat[i,]) %*% (invR_prop_phi %*% wMat[i,])
        den <- detR_trainingdata - (1/(2 * parameterMatrix[i, 3])) * t(wMat[i,]) %*% (invR_trainingdata %*% wMat[i,])
        alpha <- min(0, num - den)
      }
      
      uniff <- runif(1)
      if (alpha > log(uniff)) {
        invR_trainingdata <- invR_prop_phi
        detR_trainingdata <- det_R_prop_phi
        parameterMatrix[i, 4] <- proposed_phi
      } else {
        parameterMatrix[i, 4] <- parameterMatrix[i-1, 4]
      }
    }
    ptFinal <- proc.time() - pt
    ptFinal
    
    # acceptance rate for Beta1
    acceptRate_beta1 <- 1 - mean(parameterMatrix[-n.iteration, 1] == parameterMatrix[-1, 1])
    acceptRate_beta1
    # acceptance rate for Beta2
    acceptRate_beta2 <- 1 - mean(parameterMatrix[-n.iteration, 2] == parameterMatrix[-1, 2])
    acceptRate_beta2
    # acceptance rate for Sigma2 (Gibbs sampling/ acceptance rate should be 1)
    acceptRate_sigma2 <- 1 - mean(parameterMatrix[-n.iteration, 3] == parameterMatrix[-1, 3])
    acceptRate_sigma2
    # acceptance rate for Phi
    acceptRate_phi <- 1 - mean(parameterMatrix[-n.iteration, 4] == parameterMatrix[-1, 4])
    acceptRate_phi
    # acceptance rate for w
    acceptRate_wMat <- 1 - mean(wMat[-n.iteration, 4] == wMat[-1, 4])
    acceptRate_wMat
    
    # Prediction_MCMC
    MCMCsigma2samples <- parameterMatrix[, 3]
    MCMCphisamples <- parameterMatrix[, 4]
    exp_prediction_MCMC <- matrix(rep(NA, (0.2 * n) * n.iteration), nrow = (0.2 * n), ncol = n.iteration)
    eta_trainingdata_MCMC <- matrix(rep(NA, (0.8 * n) * n.iteration), nrow = (0.8 * n), ncol = n.iteration)
    eta_testdata_MCMC <- matrix(rep(NA, (0.2 * n) * n.iteration), nrow = (0.2 * n), ncol = n.iteration)
    CovMat_prediction_MCMC<-matrix(rep(NA,(0.2*n)^2*n.iteration), nrow=(0.2*n)^2, ncol=n.iteration)
    
    for (i in 1:n.iteration) {
      if (i %% 1000 == 0) { print(i) }
      Sigma_11 <- MCMCsigma2samples[i] * exp(-distMat[c(1:(0.8 * n)), c(1:(0.8 * n))] / MCMCphisamples[i])
      chol_inv_Sigma_11 <- chol(Sigma_11)
      inv_Sigma_11 <- chol2inv(chol_inv_Sigma_11)
      Sigma_21 <- MCMCsigma2samples[i] * exp(-distMat[c((0.8 * n + 1):n), c(1:(0.8 * n))] / MCMCphisamples[i])
      Sigma_22<-MCMCphisamples[i]*exp(-distMat[c((0.8*n+1):(n)),c((0.8*n+1):(n))]/MCMCphisamples[i])
      Sigma_12<-MCMCphisamples[i]*exp(-distMat[c(1:(0.8*n)),c((0.8*n+1):(n))]/MCMCphisamples[i])
      eta_trainingdata_MCMC[, i] <- X_trainingdata %*% parameterMatrix[i, c(1, 2)] + wMat[i,]
      eta_testdata_MCMC[, i] <- X_testdata %*% parameterMatrix[i, c(1, 2)] + Sigma_21 %*% (inv_Sigma_11 %*% wMat[i,])
      exp_prediction_MCMC[, i] <- exp(eta_testdata_MCMC[, i]) / (1 + exp(eta_testdata_MCMC[, i]))
      # Prediction 
      CovMat_prediction_MCMC[,i]<-as.numeric(Sigma_22-Sigma_21%*%inv_Sigma_11%*%Sigma_12)
    }
    
    # Final calculations
    exp_prediction_MCMC <- apply(exp_prediction_MCMC, 1, mean)
    obs_RMSPE_MCMC <- sqrt(mean((obs_testdata - exp_prediction_MCMC)^2))
    obs_RMSPE_MCMC
    rocBinary <- roc(obs_testdata, exp_prediction_MCMC)
    
    # Get the full AUC
    aucVal <- auc(rocBinary)
    aucVal
    
    burninperiod <- ceiling(n.iteration * 0.1)
    
    # Save results of the current simulation
    save(parameterMatrix, ptFinal, obs_RMSPE_MCMC, aucVal, exp_prediction_MCMC,
         acceptRate_beta1, acceptRate_beta2, acceptRate_sigma2, acceptRate_phi,
         n.iteration, burninperiod, wMat, acceptRate_wMat,
         eta_trainingdata_MCMC, eta_testdata_MCMC,CovMat_prediction_MCMC,
         file = paste0("B_phi",phiselect,"_Gen500_",sim,"_MCMC.RData"))
 

################################################################################
################################################################################

