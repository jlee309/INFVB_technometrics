################################################################################
################################################################################
################################################################################
################################################################################
# 11112024
# INLA_Count_Full
getwd()
#dev.off()
rm(list=ls())
getwd()
library(cli);library(viridis);library(fields);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm);
library(invgamma);library(pROC);library(Matrix);library(emulator);library(INLA)


################################################################################
################################################################################
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
set.seed(123)  # Set an initial seed
base_seed <- sample.int(1e6, 1)  # Randomly select a base seed


for(phiselect in 1:length(phiset) ){
  cat("Running phiselect", phiselect, "\n")
  for (sim in 1:n_simulations) {
    cat("Running simulation", sim, "\n")
    file_name <- paste0("../01DataGeneration/C_phi", phiselect, "_Gen500_", sim, ".RData")
    load(file_name)
    randomseed <- base_seed + sim + phiselect * 1000
    set.seed(randomseed)
    

################################################################################
## (1) Creating Mesh using Coordinates
mesh <- inla.mesh.2d(X_trainingdata, max.edge=c(0.1,0.2)) #(!) Let's choose the best one
plot(mesh)
dim(mesh$loc)
points(x=mesh$loc[,1], y=mesh$loc[,2],col="blue",pch=16,cex=1)
points(x=X_trainingdata[,1], y=X_trainingdata[,2],col="red",pch=16,cex=1) # Pointing out the training data

################################################################################
## (2) SPDE object
spde = inla.spde2.matern(mesh=mesh, alpha=2) #(!) Alpha check

################################################################################
## (3) the projector matrix (sparse matrix which is eta_star)
AMat <- inla.spde.make.A(mesh, loc=X_trainingdata)  
dim(AMat)
AMatCV <- inla.spde.make.A(mesh, loc=X_testdata)  
dim(AMatCV)

################################################################################
## Before (4) We might need to make an index function but not here. 
## (4) Stack all the Relevant pieces together 
## (4)-1 For the estimation 
stk.e <- inla.stack(tag='est', ## tag
                    data=list(y=obs_trainingdata), ## response
                    A=list(AMat, 1), ## two projector matrix
                    effects=list(## two elements:
                      s=1:spde$n.spde, ## RF index
                      data.frame( x1=X_trainingdata[,1],x2=X_trainingdata[,2])))

## (4)-2 Prediction fitting
stk.pred <- inla.stack(tag='pred', A=list(AMatCV, 1), 
                       data=list(y=NA), ## response as NA
                       effects=list(s=1:spde$n.spde, 
                                    data.frame(x1=X_testdata[,1],x2=X_testdata[,2])))

################################################################################
#(5) Finally Join all the stack objects p.222
stk.full <- inla.stack(stk.e, stk.pred)

################################################################################
#(6) We need to define the formula 
formula <- y ~ 0  + x1 + x2 +## fixed part
  f(s, model=spde) ## RF term

################################################################################
#(7) call INLA function 
link = c(rep(NA, dim(X_trainingdata)[1]),rep(1, dim(X_testdata)[1]))

print("Begin INLA")
res<-inla(formula,
          data=inla.stack.data(stk.full, spde=spde),
          family="poisson",  #binomial, poisson
          control.predictor=list(A=inla.stack.A(stk.full), compute=TRUE, link=link),
          control.compute = list(return.marginals.predictor=TRUE),
          keep=FALSE, verbose=TRUE)
summary(res) # Summary
print("End INLA")
print(res$cpu.used)
# save(res,
#      file=paste("ResultsStatBinary_",seqNum,".RData",sep=""), version=2)
length(res$marginals.linear.predictor) #(Q) Why is the length 332?

################################################################################
#(8) call INLA function 

field_est_mean =
  res$summary.linear.predictor[inla.stack.index(stk.full,"est")$data, "mean"]
field_est_sd =
  res$summary.linear.predictor[inla.stack.index(stk.full,"est")$data, "sd"]

# Cross Validation
field_pred_mean =
  res$summary.linear.predictor[inla.stack.index(stk.full,"pred")$data, "mean"]
field_pred_sd =
  res$summary.linear.predictor[inla.stack.index(stk.full,"pred")$data, "sd"]

lambda_prediction_INLA<-exp(field_pred_mean)#Prediction using log Link
obs_RMSPE_INLA<-sqrt(mean((obs_testdata-lambda_prediction_INLA)^2)) #CVRMSPE to be used as gold standard
lambda_RMSPE_INLA<- sqrt(mean((lambda_testdata-lambda_prediction_INLA)^2))
comptTime_INLA<-res$cpu.used # Computation Time 



################################################################################
################################################################################
#(9) Plot the linear predictor (x1*beta+w1) 
#(9)-1 Plotting Beta Distn
# Extract posterior marginals for beta
beta1_posterior_INLA <- res$marginals.fixed$x1
beta2_posterior_INLA <- res$marginals.fixed$x2
# Plot the posterior distribution for beta1
plot(beta1_posterior_INLA, type='l', main="Posterior Distribution for Beta 1", xlab="Beta 1", ylab="Density")
# Plot the posterior distribution for beta2
plot(beta2_posterior_INLA, type='l', main="Posterior Distribution for Beta 2", xlab="Beta 2", ylab="Density")
# If you want to compute summary statistics for the posterior distribution
beta1_mean <- inla.emarginal(function(x) x, beta1_posterior_INLA)
beta2_mean <- inla.emarginal(function(x) x, beta2_posterior_INLA)
cat("Posterior Mean of Beta1:", beta1_mean, "\n")
cat("Posterior Mean of Beta2:", beta2_mean, "\n")


#(9)-2 Plotting the Linear predictors
# Density of the Linear predictors
xSeq<-seq(-3,3,length.out=1000)
ySeq<-dnorm(x = xSeq , mean = field_pred_mean[1] , sd = field_pred_sd[1])
plot(x=xSeq , y = ySeq , typ="l" , main = "Density for Linear Predictor 1")


#(9)-3 Spatial Random Effects W (But not sure this is correct way or not)
# For the training data
# Extract the index for the spatial random effects for the training data
spatial_index_trainingdata <- inla.stack.index(stk.full, "est")$data
# Extract posterior summaries for the spatial random effects
w_trainingdata_INLA_mean <- res$summary.random$s$mean[spatial_index_trainingdata]
w_trainingdata_INLA_SD <- res$summary.random$s$sd[spatial_index_trainingdata]


# For the test data
# Extract the index for the spatial random effects for the test data
spatial_index_testdata<-inla.stack.index(stk.full, "pred")$data
# Extract posterior summaries for the spatial random effects
w_testdata_INLA_mean <- res$summary.random$s$mean[spatial_index_testdata]
w_testdata_INLA_SD <- res$summary.random$s$sd[spatial_index_testdata]


################################################################################
#(10) Save 
save( obs_RMSPE_INLA ,comptTime_INLA,lambda_prediction_INLA, lambda_RMSPE_INLA,
     beta1_posterior_INLA, beta2_posterior_INLA,
     w_trainingdata_INLA_mean,w_trainingdata_INLA_SD,
     w_testdata_INLA_mean,w_testdata_INLA_SD,
     field_pred_mean,field_pred_sd,
     field_est_mean,field_est_sd,
     file = paste0("C_phi",phiselect,"_Gen500_",sim,"_INLA.RData"))
  }
}
