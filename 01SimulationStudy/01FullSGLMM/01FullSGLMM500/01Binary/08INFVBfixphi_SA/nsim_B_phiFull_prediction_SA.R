################################################################################
################################################################################
################################################################################
################################################################################
# 10142024
# nsim_SA_INFVBPhi
rm(list=ls())
library(mvtnorm);
library(cli);library(viridis);library(fields);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm)
library(invgamma);library(pROC);library(Matrix);library(emulator)

################################################################################
################################################################################
# Sampling Method for estimating Eta and predicting Eta 
# By using Multivariate Normal for v(Beta, W) 

################################################################################
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name1 <- paste0("../01DataGeneration/B_phi", phiselect, "_Gen500_", sim, ".RData")
    load(file_name1)
    file_name2 <- paste0("../03INFVBfixphi_ParallelComputing/B_phi", phiselect, "_Gen500_", sim, "_VBphiParallel.RData")
    load(file_name2)
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)


################################################################################
# (1) Getting 50 samples for V(Beta, W) and ssq
numofsamples<-50
sample_V_list <- list()
sample_ssq_list <- list()

## Multivariate normal distribution for w 
for (k in 1:iter.thetaphi) {
  
  if(k%%100 == 0) {print(k)}
  #Setting the covariance matrix for v(beta and W)
  thetaC_covMatV_Mat<-diag(thetaC_covMatV[k,])
  
  # Draw samples
  sample_V<-rmvnorm(numofsamples, mean=thetaC_MeanV[k,], sigma= thetaC_covMatV_Mat)
  sample_ssq <- rinvgamma(numofsamples, shape = thetaC_sigma2[k, 1], rate = thetaC_sigma2[k, 2])
  
  
  # Save samples into lists
  sample_V_list[[k]] <- sample_V
  sample_ssq_list[[k]] <- sample_ssq
}

#Just for Checking 
dim(sample_V_list[[1]])
sample_ssq_list[[1]]

################################################################################
# (2) Getting numofsamples Eta for each discretization
################################################################################
sample_eta_trainingdata_VB <- list()
sample_eta_testdata_VB <- list()
sample_CovMat_prediction_VB <- list()
sample_W_test_list<-list()

for(i in 1:iter.thetaphi){
  if(i%%100 == 0) {print(i)}
  
  # Initialize these lists at the i-th position before using them
  sample_eta_trainingdata_VB[[i]] <- list()  # Initialize an empty list for the i-th position
  sample_eta_testdata_VB[[i]] <- list()      # Same for test data
  sample_CovMat_prediction_VB[[i]] <- list() # Same for prediction covariance matrix
  
  for(j in 1:numofsamples){
    Sigma_11 <- sample_ssq_list[[i]][j] * exp(-distMat[c(1:(0.8 * n)), c(1:(0.8 * n))] / thetaphi[i])
    inv_Sigma_11 <- solve(Sigma_11)
    Sigma_21 <- sample_ssq_list[[i]][j] * exp(-distMat[c((0.8 * n + 1):n), c(1:(0.8 * n))] / thetaphi[i])
    Sigma_22 <- sample_ssq_list[[i]][j] * exp(-distMat[c((0.8 * n + 1):n), c((0.8 * n + 1):n)] / thetaphi[i])
    Sigma_12 <- sample_ssq_list[[i]][j] * exp(-distMat[c(1:(0.8 * n)), c((0.8 * n + 1):n)] / thetaphi[i])
    
    sample_eta_trainingdata_VB[[i]][[j]] <- X_trainingdata %*%sample_V_list[[i]][j,][c(1:p)]+sample_V_list[[i]][j,][c((p+1):(0.8 * n+p))]
    sample_CovMat_prediction_VB_ing<-Sigma_22 - Sigma_21 %*% inv_Sigma_11 %*% Sigma_12
    sample_CovMat_prediction_VB[[i]][[j]] <- as.numeric(sample_CovMat_prediction_VB_ing)
    sample_W_test_list[[i]]<-rmvnorm(n=100, mean=Sigma_21%*%(inv_Sigma_11%*%sample_V_list[[i]][j,][c((p+1):(0.8 * n+p))]), sigma =sample_CovMat_prediction_VB_ing)
    sample_eta_testdata_VB[[i]][[j]] <- X_testdata %*%sample_V_list[[i]][j,][c(1:p)]+sample_W_test_list[[i]][j,]
  }
}



# Just for checking 
dim(sample_W_test_list[[1]])
sample_eta_trainingdata_VB[[1]][1]
sample_eta_testdata_VB[[1]][1]
sample_CovMat_prediction_VB[[1]][1]


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-Delete for the memory#-#-#-#-#-#-#-#-#-#-#-#-#-#
# rm(sample_V_list, sample_ssq_list,sample_CovMat_prediction_VB,sample_W_test_list)
################################################################################
# (3) Get the ELBO*Samples_Estimation (note : Not a prediction now)
################################################################################

# Initialize an empty list to store the vectors for each i (i1, i2, etc.)
ELBOmean_list_est <- list()

# Main nested loop
for (i in 1:iter.thetaphi) {
  if(i%%100 == 0) {print(i)}
  # Initialize an empty vector to store the values for each i
  ELBOmean_vector <- c()
  
  # Loop over the k values 
  for (k in 1:(0.8 * n)) {
    mean_value <- c()
    
    # Calculate mean_value for each j
    for (j in 1:numofsamples) {
      mean_value_cal <- unlist(sample_eta_trainingdata_VB[[i]][j])[k]
      mean_value <- c(mean_value, mean_value_cal)
    }
    
    # Calculate ELBOmean_value and store it dynamically
    assign(paste0("ELBOmean_value_i", i, "_k", k), mean(mean_value) * ELBOvectorweights[i])
    
    # Retrieve the dynamically assigned scalar value and append to the vector
    dynamic_name <- paste0("ELBOmean_value_i", i, "_k", k)
    scalar_value <- get(dynamic_name)
    ELBOmean_vector <- c(ELBOmean_vector, scalar_value)
    
    # Optionally remove the dynamically created object after use
    rm(list = dynamic_name)
  }
  
  # Store the vector in the list
  ELBOmean_list_est[[paste0("ELBOmean_vector_i", i)]] <- ELBOmean_vector
}



#just for checking 
ELBOmean_list_est[[1]]
length(ELBOmean_list_est[[1]])

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-Delete for the memory#-#-#-#-#-#-#-#-#-#-#-#-#-#
# rm(sample_eta_trainingdata_VB)

################################################################################
# (4) Making the sample covariance matrix for the estimation 
################################################################################

cov_sample_eta_trainingdata_VB_all<-matrix(data=rep(0,(0.8*n)^2), nrow=(0.8*n), ncol=(0.8*n))
for(i in 1:iter.thetaphi){
  cov_sample_eta_trainingdata_VB<-c()
  for(j in 1:numofsamples){
    cov_sample_eta_trainingdata_VB<-
      cbind(cov_sample_eta_trainingdata_VB, unlist(sample_eta_trainingdata_VB[[i]][j]))
  }
  cov_sample_eta_trainingdata_VB<-t(cov_sample_eta_trainingdata_VB)
  cov_sample_eta_trainingdata_VB<-cov(cov_sample_eta_trainingdata_VB)
  assign(paste0("cov_sample_eta_trainingdata_VB_i",i),cov_sample_eta_trainingdata_VB*ELBOvectorweights[i] )
  dynamic_name<-paste0("cov_sample_eta_trainingdata_VB_i",i)
  cov_sample_eta_trainingdata_VB_i<-get(dynamic_name)
  cov_sample_eta_trainingdata_VB_all<-cov_sample_eta_trainingdata_VB_all+cov_sample_eta_trainingdata_VB_i
  # Optionally remove the dynamically created object after use
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-Delete for the memory#-#-#-#-#-#-#-#-#-#-#-#-#-#
  rm(list = dynamic_name)
}

#Just for checking
dim(cov_sample_eta_trainingdata_VB_all)

################################################################################
# (5) Plotting Etas in the training datasets
################################################################################
# The real Data
eta_trainingdata_real <-X_trainingdata%*%beta+w_trainingdata
eta_testdata_real <-X_testdata%*%beta+w_testdata
################################################################################
# Making mu_sample_est_VB
################################################################################
mu_sample_est_VB<-Reduce("+",ELBOmean_list_est)
mu_sample_est_VB

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-Delete for the memory#-#-#-#-#-#-#-#-#-#-#-#-#-#
# rm(ELBOmean_list_est)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Prediction 
################################################################################
# (3-1) Get the ELBO*Samples_Prediction
################################################################################

# Initialize an empty list to store the vectors for each i (i1, i2, etc.)
ELBOmean_list_pred <- list()

# Main nested loop
for (i in 1:iter.thetaphi) {
  if(i%%100 == 0) {print(i)}
  # Initialize an empty vector to store the values for each i
  ELBOmean_vector <- c()
  
  # Loop over the k values (1 to 20)
  for (k in 1:(0.2 * n)) {
    mean_value <- c()
    
    # Calculate mean_value for each j
    for (j in 1:numofsamples) {
      mean_value_cal <- unlist(sample_eta_testdata_VB[[i]][j])[k]
      mean_value <- c(mean_value, mean_value_cal)
    }
    
    # Calculate ELBOmean_value and store it dynamically
    assign(paste0("ELBOmean_value_i", i, "_k", k), mean(mean_value) * ELBOvectorweights[i])
    
    # Retrieve the dynamically assigned scalar value and append to the vector
    dynamic_name <- paste0("ELBOmean_value_i", i, "_k", k)
    scalar_value <- get(dynamic_name)
    ELBOmean_vector <- c(ELBOmean_vector, scalar_value)
    
    # Optionally remove the dynamically created object after use
    rm(list = dynamic_name)
  }
  
  # Store the vector in the list
  ELBOmean_list_pred[[paste0("ELBOmean_vector_i", i)]] <- ELBOmean_vector
}

#just for checking
ELBOmean_list_pred[[1]]

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-Delete for the memory#-#-#-#-#-#-#-#-#-#-#-#-#-#
# rm(sample_eta_testdata_VB)
################################################################################
# (4-1) Making the sample covariance matrix for the estimation 
################################################################################

cov_sample_eta_testdata_VB_all<-matrix(data=rep(0,(0.2*n)^2), nrow=(0.2*n), ncol=(0.2*n))
for(i in 1:iter.thetaphi){
  cov_sample_eta_testdata_VB<-c()
  for(j in 1:numofsamples){
    cov_sample_eta_testdata_VB<-
      cbind(cov_sample_eta_testdata_VB, unlist(sample_eta_testdata_VB[[i]][j]))
  }
  cov_sample_eta_testdata_VB<-t(cov_sample_eta_testdata_VB)
  cov_sample_eta_testdata_VB<-cov(cov_sample_eta_testdata_VB)
  assign(paste0("cov_sample_eta_testdata_VB_i",i),cov_sample_eta_testdata_VB*ELBOvectorweights[i] )
  dynamic_name<-paste0("cov_sample_eta_testdata_VB_i",i)
  cov_sample_eta_testdata_VB_i<-get(dynamic_name)
  cov_sample_eta_testdata_VB_all<-cov_sample_eta_testdata_VB_all+cov_sample_eta_testdata_VB_i
  # Optionally remove the dynamically created object after use
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-Delete for the memory#-#-#-#-#-#-#-#-#-#-#-#-#-#
  rm(list = dynamic_name)
}

#Just for checking
dim(cov_sample_eta_testdata_VB_all)

################################################################################
# Making mu_sample_est_VB
################################################################################
mu_sample_pred_VB<-Reduce("+",ELBOmean_list_pred)
mu_sample_pred_VB
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-Delete for the memory#-#-#-#-#-#-#-#-#-#-#-#-#-#
# rm(ELBOmean_list_pred)

################################################################################
################################################################################
save(eta_trainingdata_real, eta_testdata_real,
     cov_sample_eta_trainingdata_VB_all, cov_sample_eta_testdata_VB_all,
     mu_sample_est_VB,mu_sample_pred_VB,
     file = paste0("SA_B_phi",phiselect,"_Gen500_",sim,"_Full_INFVBphi.RData"))
  }
}