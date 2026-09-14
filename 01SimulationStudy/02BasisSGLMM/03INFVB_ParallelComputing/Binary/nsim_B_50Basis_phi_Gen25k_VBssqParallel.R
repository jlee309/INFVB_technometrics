################################################################################
################################################################################
# INFVB for Binary Data using 50 Basis functions 
# when sigma2 is discretized (Parallelization Computing)  
# Rcpp functions are used 
# Only using Diag for the covariance matrix of V
# 10242024 
################################################################################
################################################################################
rm(list=ls())
getwd()
################################################################################
################################################################################
library(inline);library(cli);library(viridis);library(fields);library(pgdraw);
library(MASS);library(tictoc);library(invgamma);library(mvtnorm);
library(pROC);library(Matrix);library(emulator);library(parallel);
library(Rcpp);library(RcppEigen) 

# sourceCpp("matsourceEigen.cpp")
sourceCpp("matsourceEigen.cpp", cacheDir = "/scratch/jlee309/rcpp-cache/Basis_INFVB_phi")
detectCores() 
nprocs <-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
this_cluster <- parallel::makeForkCluster(nprocs) 

################################################################################
################################################################################
## Let the cache directory work sequentially.  
# for (i in seq_along(this_cluster)) {
#   clusterEvalQ(this_cluster[i], {
#     sourceCpp("matsourceEigen.cpp", cacheDir = "/scratch/jlee309/rcpp-cache/Basis_INFVB_phi")
#   })
# }
################################################################################
################################################################################
# Fast Self Made Functions
quad.diag.New<-function (M, x) { # Diagonal elements of a quadratic form diag(t(x)%*%M%*%x)
  colSums(crossprod(M, x) * x)
}
quad.form.new<-function(M,x){ # t(x)%*%M%*%x
  crossprod(crossprod(M, x), x)
}
computeTrace<-function(A,B){  #sum(diag(A%*B))
  sum(colSums(t(A)*B))
}
################################################################################
################################################################################
# Number of simulations
n_simulations <- 50 #Number of simulations per each scenario
phiset=c(0.1, 0.3, 0.5, 0.7)
# phiset<-1
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
    
    
    # Set Basis Functions
    ################################################################################
    ################################################################################
    numofbasis<-50 
    M<-refEigen_trainingdata[,1:numofbasis]
    M_CV<-refEigen_testdata[,1:numofbasis]
    ################################################################################
    ################################################################################
    # Variational Bayes (VB)
    ## Preliminaries
    ## Initial Settings for VB
    p<-length(beta)
    n<-length(obsBin_trainingdata)
    iter<-1000
    
    
    ## Beta
    VB_mean_beta<-c(0,0)
    VB_covMat_beta<-100*diag(2)
    VB_inv_Sigma_beta_prior<-solve(VB_covMat_beta)
    
    
    ## Prior Distribution For Sigma2 (Inverse Gamma) 
    VB_alpha_sigma<-0.1 ; VB_beta_sigma<-0.1
    #VB_E_sigma2Inv<-VB_iter[1,1]/VB_iter[1,2]
    
    
    ## Prior Distribution For Delta 
    VB_mu_delta<-c(rep(0,numofbasis)) ;VB_sigma_delta<-1*diag(numofbasis)
    VB_cholD<-chol(VB_sigma_delta)
    VB_inv_sigma_delta<-chol2inv(VB_cholD)
    
    
    ## Making C matrix 
    dim(refEigen_trainingdata)
    C<-cbind(X_trainingdata,M)
    tC<-t(C)
    dim(C)
    
    #Making V(beta+W) Matrix
    VB_MeanbetaMat<-matrix(NA, nrow=iter, ncol=p) 
    dim(VB_MeanbetaMat)
    VB_mean_deltaMat<-matrix(NA,nrow=iter,ncol=(numofbasis)) 
    dim(VB_mean_deltaMat)
    
    
    #First Line
    VB_MeanbetaMat[1,]<-c(rep(1,p))
    VB_mean_deltaMat[1,]<-c(rep(0,numofbasis))
    
    
    #Making VB_MeanV Matrix
    VB_MeanV<-cbind(VB_MeanbetaMat, VB_mean_deltaMat)
    dim(VB_MeanV)
    
    
    #Making covMatV Matrix 
    VB_covMatV<-matrix(NA, nrow=iter, ncol=(numofbasis)+p)
    dim(VB_covMatV)
    VB_covMatV[1,]<-rep(1,((numofbasis)+p))
    zeroMatpn<-matrix(0, nrow=p, ncol=(numofbasis)) 
    zeroMatnp<-matrix(0, nrow=(numofbasis), ncol=p)
    
    
    #Prior
    invVMatmaking1<-cbind(VB_inv_Sigma_beta_prior,zeroMatpn)
    dim(invVMatmaking1)
    invVMatmaking2<-cbind(zeroMatnp,1*VB_inv_sigma_delta)
    dim(invVMatmaking2)
    invVMat<-rbind(invVMatmaking1,invVMatmaking2)
    dim(invVMat)
    
    
    # Observations
    obsMinHalf<-obsBin_trainingdata-0.5
    tCobsMinHalf<-as.numeric(tC%*%obsMinHalf)
    
    ################################################################################
    ################################################################################
    #Making temporary lambda to get xi (Using Rcpp Functions)
    #Method A-1(Without Using Rcpp)
    #system.time({xi<-sqrt(quad.diag.New(M=(VB_covMatV[1,]+VB_MeanV[1,]%*%t(VB_MeanV[1,])),x = tC))})
    
    #Method A-2(With Using Rcpp) #Method A-2 is much faster than Method A-1
    system.time({
      foo<-(VB_covMatV[1,]+VB_MeanV[1,]%*%t(VB_MeanV[1,]))
      xi<-(sqrt(myfun3(foo,tC))) 
    })
    
    lambda<- -tanh(xi/2)/(4*xi)
    lambda<-Diagonal(x = lambda)
    
    #Method B-1(Without Using Rcpp) 
    #system.time({tCLambdaC_1<-quad.form.new(M=lambda,x=C)}) 
    
    #Method B-2(With Using Rcpp) #Method B-2 is much faster than B-1
    system.time({
      baz<-as.matrix(lambda%*%C)
      tCLambdaC<-myfun2(C,baz)}) 
    
    
    ################################################################################
    ################################################################################
    #INFVB Preparation
    
    ## ELBOvector 
    ELBOvector<-ELBOjthcal<-numeric()
    ELBOjthcal[1]<-10
    
    #INFVB Discretizing ThetaD:{sigma2} 
    iter.thetaphi=1000 #Number of Discretizing ThetaD 
    x_sigma2=seq(0.001,2000,length.out=iter.thetaphi)
    theta_sigma2=dinvgamma(x=x_sigma2,shape=VB_alpha_sigma, rate=VB_beta_sigma)
    #plot(x_sigma2,theta_sigma2,col="blue")
    
    ## Final Saving Matrix for V(Beta, W) and sigma2
    thetaC_MeanV<-matrix(NA, nrow=iter.thetaphi, ncol=(numofbasis+p))
    dim(thetaC_MeanV)
    thetaC_covMatV<-matrix(NA, nrow=iter.thetaphi, ncol=(numofbasis+p))
    dim(thetaC_covMatV)
    
    
    VB_covMatchol<-chol(invVMat -2*tCLambdaC)
    VB_covMatVing<-chol2inv(VB_covMatchol)
    VB_covMatVing<-as.matrix(VB_covMatVing)
    VB_covMatV[1,]<- diag(VB_covMatVing)
    VB_MeanV[1,] <- as.numeric(VB_covMatVing%*%tCobsMinHalf)
    
    
    m<-numofbasis
    basisInd<-(1:m)+p
    betaInd<-1:p
    
    
    ################################################################################ 
    ################################################################################
    #VB Start
    VB_pt<-proc.time()
    
    parallel::clusterExport( this_cluster, varlist=c("zeroMatnp",
                                                     "x_sigma2",
                                                     "m",
                                                     "invVMatmaking1",
                                                     "theta_sigma2",
                                                     "tCLambdaC",
                                                     "tCobsMinHalf",
                                                     "tC","C",
                                                     "quad.form.new", 
                                                     "ELBOjthcal",
                                                     "VB_covMatV",
                                                     "VB_MeanV"
    ))
    
    outputMat_list<-parallel::clusterApply(this_cluster, 1:iter.thetaphi, function (j) {
      
      # Updatng invVMat 
      invVMatmaking2<-cbind(zeroMatnp,Diagonal(x=rep((1/x_sigma2[j]),m)))
      invVMat<-rbind(invVMatmaking1,invVMatmaking2)
      invVMat_chol<-chol(invVMat)
      # VMat<-chol2inv(invVMat_chol)
      # VMat_chol<-chol(VMat)
      
      ELBO3<-log(theta_sigma2[j])
      ELBOjthcal[1]<-10
      k=1
      bar<-77 #Arbitrary number 
      
      while(bar>exp(1e-3)) {
        k<-k+1
        print(k)
        # Still Updatng invVMat 
        VB_covMatchol<-chol(invVMat -2*tCLambdaC)
        VB_covMatVing<-chol2inv(VB_covMatchol)
        VB_covMatVing<-as.matrix(VB_covMatVing)
        VB_covMatV[k,]<- diag(VB_covMatVing)
        VB_MeanV[k,] <- as.numeric(VB_covMatVing%*%tCobsMinHalf)
        
        #Update xi- Computational Bottlenecks with Matrix -Matrix quadratic forms
        foo<-(VB_covMatVing+VB_MeanV[k-1,]%*%t(VB_MeanV[k-1,]))
        xi<-sqrt(myfun3(foo,tC)) 
        lambda<- -tanh(xi/2)/(4*xi)
        sum_psi_xi<-sum(xi*0.5-log(1+exp(xi))+xi*(tanh(xi/2)/(4)))
        lambda<-Diagonal(x = lambda)
        #lambda<-diag(x = lambda)
        baz<-as.matrix(lambda%*%C)
        tCLambdaC<-myfun2(C,baz)
        
        
        ################################################################################
        # Stopping Criteria (By using difference between ELBO)
        #ELBO preparation
        ELBO1<-sum(quad.form.new(M =tCLambdaC , x= VB_MeanV[k,]),
                   sum(diag(tCLambdaC%*%VB_covMatVing)),
                   as.numeric(-0.5*quad.form.new(M = invVMat, x= VB_MeanV[k,])),
                   -0.5*sum(diag(invVMat%*%VB_covMatVing)), 
                   sum_psi_xi)
        
        ELBO2<- tCobsMinHalf%*%VB_MeanV[k,]
        ELBO4<- -sum(log(diag(VB_covMatchol))) #same as 1/2*log(det(VB_covMatVing))
        ELBO5<- sum(log(diag(invVMat_chol))) #same as -sum(log(diag(VMat_chol)))
        ELBOjthcal[k] <- ELBO1+ELBO2+ELBO3+ELBO4+ELBO5
        
        #Making ELBO for the Stopping Criteria
        bar<-abs(ELBOjthcal[k]-ELBOjthcal[k-1])
      }
      return(c(ELBOjthcal[k],VB_MeanV[k,],VB_covMatV[k,]))
    })
    VB_ptFinal<-proc.time()-VB_pt
    VB_ptFinal
    # parallel::stopCluster(cl)
    
    
    outputMat<-matrix(unlist(outputMat_list), ncol=length(outputMat_list))
    dim(outputMat)
    ELBOvector<-outputMat[1,]
    
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    # Exponential ELBO
    weights<-exp(ELBOvector-max(ELBOvector))
    ELBOvectorweights<-weights/sum(weights)
    
    ################################################################################
    ################################################################################
    
    for ( i in 1:iter.thetaphi){
      thetaC_MeanV[i,]<-outputMat[2:((numofbasis)+p+1),i]
    }
    
    for ( i in 1:iter.thetaphi){
      thetaC_covMatV[i,]<-outputMat[((numofbasis)+p+2):(1+p+numofbasis+(p+numofbasis)),i]
    }
    
    
    ############################################################################
    ############################################################################
    # Making new Density according to ELBOvector rates
    ############################################################################
    ############################################################################
    # Beta1 Weighted Density
    par(mfrow=c(1,1))
    matDensityBeta1<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
    xSeqb1<-seq(-4,4, length.out=iter.thetaphi)
    for ( k in 1:iter.thetaphi){
      densitiesBeta1 <-dnorm(x=xSeqb1, mean=thetaC_MeanV[k,1], sd=sqrt(thetaC_covMatV[k,1]))
      matDensityBeta1[,k]<-densitiesBeta1
    }
    dim(matDensityBeta1)
    newdensityBeta1<-matDensityBeta1%*%ELBOvectorweights
    #plot(x=xSeqb1, y=newdensityBeta1, typ="l", pch=16, col="red",main="Beta1")
    
    # Beta2 Weighted Density
    matDensityBeta2<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
    xSeqb2<-seq(-4,4, length.out=iter.thetaphi)
    for ( k in 1:iter.thetaphi){
      densitiesBeta2 <-dnorm(x=xSeqb2, mean=thetaC_MeanV[k,p], sd=sqrt(thetaC_covMatV[k,2]))
      matDensityBeta2[,k]<-densitiesBeta2
    }
    dim(matDensityBeta2)
    newdensityBeta2<-matDensityBeta2%*%ELBOvectorweights
    dim(newdensityBeta2)
    #plot(x=xSeqb2, y=newdensityBeta2, typ="l", pch=16, col="red", xlim=c(0,2), main="Beta2")
    
    
    ################################################################################
    ################################################################################
    ## delta
    ################################################################################
    ################################################################################
    # delta3 Weighted Density
    q=3
    matDensitydelta3<-matrix(NA,ncol=iter.thetaphi,nrow=iter.thetaphi)
    xSeqdelta3<-seq(-20,20, length.out=iter.thetaphi)
    for ( k in 1:iter.thetaphi){
      densitiesdelta3 <-dnorm(x=xSeqdelta3, mean=thetaC_MeanV[k,(p+q)], sd=sqrt(thetaC_covMatV[k,p+q]))
      matDensitydelta3[,k]<-densitiesdelta3
    }
    dim(matDensitydelta3)
    newdensitydelta3<-matDensitydelta3%*%ELBOvectorweights
    
    
    ################################################################################
    ################################################################################
    ## Cross Validation
    ################################################################################
    ################################################################################
    #Prediction_VB
    ################################################################################
    n<-dim(M)[1]+dim(M_CV)[1]
    exp_prediction_VB<-matrix(rep(NA,(0.2*n)*iter.thetaphi), nrow=(0.2*n), ncol=iter.thetaphi)
    eta_testdata_VB<-matrix(rep(NA,(0.2*n)*iter.thetaphi), nrow=(0.2*n), ncol=iter.thetaphi)
    eta_trainingdata_VB<-matrix(rep(NA,(0.8*n)*iter.thetaphi), nrow=(0.8*n), ncol=iter.thetaphi)
    
    dim(exp_prediction_VB)
    for(i in 1:iter.thetaphi){
      if(i%%100==0) {print(i)}
      eta_trainingdata_VB[,i]<-X_trainingdata%*%thetaC_MeanV[i,c(1,2)]+M%*%thetaC_MeanV[i,c((1+p):(p+numofbasis))]
      eta_testdata_VB[,i]<-X_testdata%*%thetaC_MeanV[i,c(1,2)]+M_CV%*%thetaC_MeanV[i,c((1+p):(p+numofbasis))]
      exp_prediction_VB[,i]<- exp(eta_testdata_VB[,i])/(1+exp(eta_testdata_VB[,i]))
    }
    dim(exp_prediction_VB)
    exp_prediction_VB<-exp_prediction_VB%*%ELBOvectorweights
    exp_RMSPE_VB<-sqrt(mean((obsBin_testdata-exp_prediction_VB)^2))
    exp_RMSPE_VB
    
    library(pROC)
    rocBinary <- roc(obsBin_testdata, exp_prediction_VB)
    aucVal_VB<-auc(rocBinary)
    save(ELBOvectorweights,newdensityBeta1, newdensityBeta2,iter.thetaphi,
         newdensitydelta3,exp_RMSPE_VB,exp_prediction_VB,VB_ptFinal,x_sigma2,aucVal_VB,
         xSeqb1, xSeqb2,  xSeqdelta3,
         thetaC_MeanV, thetaC_covMatV, eta_testdata_VB, eta_trainingdata_VB,
         numofbasis,p,
         file=paste0("50Basis_B_phi",phiselect,"_Gen25k_",sim,"_VBssqParallel.RData")) 
  }
}