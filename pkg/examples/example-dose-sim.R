# This example demonstrates the simulation of a bolus 
# in a 2-compartment system.

rm(list=ls())

k1 = 0.05; k2 = 0.05; ke = 0.07;

Model.SimDose <- list(
                      Matrices = function(phi) {
                        kei <- phi[["kei"]]
                        matA <- matrix( c(-(k1+kei) ,  k2 ,   
                                                 k1 , -k2 ) ,nrow=2,byrow=T)
                        matC <- matrix(c(1,0),nrow=1)
                        list(matA=matA,matC=matC)
                      },
                      X0 = function(Time=Na,phi,U=Na) {
                        matrix(0,nrow=2)
                      },
                      SIG = function(phi) { 
                        diag( c(0,0) )
                      },
                      S = function(phi) {
                        matrix(phi[["S"]])
                      },
                      h = function(eta,theta,covar) {
                        phi <- theta
                        phi[["kei"]] <- theta[["ke"]]*exp(eta[1])
                        phi                      },
                      ModelPar = function(THETA){
                        list(theta=list(ke=THETA[1], S=THETA[2]),
                             OMEGA=matrix(THETA[3]) )
                      },
                      Dose = list(
                        Time = c(30,180),
                        State = c(1, 1),
                        Amount = c(1500,1500)
                        )
                      )

# Create Simulation Timeline 
SimDose.Subj <- 2
SimDose.Data <- vector(mode="list",length=SimDose.Subj)
for (i in 1:SimDose.Subj) 
  SimDose.Data[[i]]$Time <- seq(from=15,by=15,length=30)




#############
# Simulation
#############

detach(package:PSM)
library(PSM,lib.loc="~/PSM/Rpackages/gridterm")

#                   ke   S  OMEGA
SimDose.THETA <-  c(0.03 , 10 , 1)

SimDose.Data <- PSM.simulate(Model.SimDose, SimDose.Data, SimDose.THETA, dt=.1 , individuals=SimDose.Subj)

# View the Simulated datastructure
names(SimDose.Data[[1]])

# Plot of the simulations
par(mfcol=c(3,SimDose.Subj))
for(id in 1:SimDose.Subj) {
  plot(SimDose.Data[[id]]$Time , SimDose.Data[[id]]$Y,
         ylab="Observations", xlab=paste('individual',id))
  for(i in 1:2) {
    plot(SimDose.Data[[id]]$Time , SimDose.Data[[id]]$X[i,],type="l",
         ylab=paste('state',i), xlab=paste('individual',id,', eta:',round(SimDose.Data[[id]]$eta,3)))
    rug(SimDose.Data[[id]]$Time)
  }
}



###########
# Smoothing
###########
#source("~/PSM/PSM/R/PSM.smooth.R")
out <- PSM.smooth(Model = Model.SimDose, Data = SimDose.Data, THETA = SimDose.THETA, subsample = 10, etaList=matrix(c(0,0),nrow=1))

# View the data structure
names(out[[1]])

#Plot of the smoothed estimates
par(mfcol=c(3,SimDose.Subj))
for(id in 1:SimDose.Subj) {
  plot(SimDose.Data[[id]]$Time , SimDose.Data[[id]]$Y,
         ylab="Observations", xlab=paste('individual',id))
  for(i in 1:2) {
    plot(out[[id]]$Time , out[[id]]$Xs[i,],type="l",
         ylab=paste('smooth state',i), xlab=paste('individual',id))
    rug(SimDose.Data[[id]]$Time)
  }
}
