################################################################################
################################################################################
# 11102024 
# INLA_Basis_Binary Data
################################################################################
################################################################################
library(INLA)
# dev.off()
# Load the data
rm(list=ls())
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed

for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name1 <- paste0("../../01DataGeneration/phi",phiselect,"_Gen25k_",sim,"EigenBasis.RData")
    file_name2 <-paste0("../../01DataGeneration/phi",phiselect,"_Gen25k_",sim,"SpatialData.RData")
    load(file_name1)
    load(file_name2)
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
################################################################################
# Prepare basis functions and data
numofbasis <- 100 
M <- refEigen_trainingdata[, 1:numofbasis]
M_CV <- refEigen_testdata[, 1:numofbasis]
dim(M); dim(M_CV)
# Combine predictors and basis functions
group_trainingdata <- 1:nrow(refEigen_trainingdata)  # Ensure this covers the correct range
group_testdata <- (1:nrow(refEigen_testdata)) + max(group_trainingdata)  # Ensure no overlap
group_trainingdata_all<-c(group_trainingdata,group_testdata)

H_tilde_trainingdata <- cbind(X_trainingdata, M)
dim(H_tilde_trainingdata)
H_tilde_testdata <- cbind(X_testdata, M_CV)
dim(H_tilde_testdata)
H_tilde_all<-rbind(H_tilde_trainingdata, H_tilde_testdata)
dim(H_tilde_all)

# Set column names for consistency
colnames(H_tilde_trainingdata) <- paste0("V", 1:ncol(H_tilde_trainingdata))
colnames(H_tilde_testdata) <- paste0("V", 1:ncol(H_tilde_testdata))
colnames(H_tilde_all) <- paste0("V", 1:ncol(H_tilde_all))
yy=c(obsBin_trainingdata, rep(NA,(0.2*n)))

# Combine into a data frame for INLA
data <- data.frame(y = yy, H_tilde_all, group = group_trainingdata_all)

################################################################################
# Fit the INLA model
formula <- as.formula(paste("y ~ -1 + ", paste(colnames(H_tilde_all), collapse = " + "), " + f(group, model = 'iid')"))

model <- inla(formula, 
              family = "binomial", 
              data = data, 
              control.predictor = list(link = 1))

# Print summary of the model
summary(model)
# Print summary of the model
summary(model)
TotTime_INLA<-model$cpu.used[4]
################################################################################
# Extract estimated coefficients and linear predictors
INLA_Basis_Mean <- model$summary.fixed[, 1]
INLA_Basis_SD <- model$summary.fixed[, 2]
linear_predictor_mean <- model$summary.linear.predictor$mean
linear_predictor_sd <- model$summary.linear.predictor$sd
head(linear_predictor_mean)
head(linear_predictor_sd)
################################################################################
# Extract predicted values and actual observed values
predicted_values <- model$summary.fitted.values$mean  # Extract predicted means
actual_values <- data$y  # Actual values (including training and test data)
# Identify which values correspond to the test data (where 'y' is NA)
test_indices <- which(is.na(yy))
# Compute RMSPE for test data
predicted_test <- predicted_values[test_indices]
actual_test <- obsPois_testdata  # Replace with the correct test data vector
# RMSPE calculation
obs_RMSPE_INLA <- sqrt(mean((predicted_test - actual_test)^2, na.rm = TRUE))

library(pROC)
rocBinary <- roc(obsBin_testdata, predicted_test)
# Get the full AUC
aucVal_INLA<-auc(rocBinary)
aucVal_INLA

save(INLA_Basis_Mean, INLA_Basis_SD, 
     linear_predictor_mean,linear_predictor_sd,
     TotTime_INLA, obs_RMSPE_INLA, aucVal_INLA,
     file=paste0("100Basis_B_phi",phiselect,"_Gen25k_",sim,"_INLA.RData")) 
  }
}
