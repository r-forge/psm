# This example demonstrates the simulation of a bolus dose
# in a 2-compartment system.
V1 <- 5      #[L]
V2 <- 10     #[L]
CLd <- 0.005 #[L/min]
CL <- 0.002  #[L/min]

Model.SimDose <- list(
                      Matrices = function(phi) {
                        V1i <- phi[["V1i"]]
                        matA <- matrix(c(-(CL+CLd)*V1i ,  CLd*V2 ,   
                                           CLd*V1i , -CLd*V2 ) ,nrow=2,byrow=T)
                        matC <- matrix(c(1/V1i,0),nrow=1)
                        list(matA=matA,matC=matC)
                      },
                      X0 = function(Time=Na,phi,U=Na) {
                        matrix(0,nrow=2)
                      },
                      SIG = function(phi) { 
                        sig11 <- phi[["sig11"]]
                        matrix(c( sig11,0,
                                 -sig11,0), nrow=2, byrow=T)
                      },
                      S = function(phi) {
                        matrix(phi[["S"]])
                      },
                      h = function(eta,theta,covar) {
                        phi <- theta
                        phi[["V1i"]] <- V1*exp(eta[1])
                        phi                      },
                      ModelPar = function(THETA){
                        list(theta=list(sig11=THETA[1], S=THETA[2]),
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
  SimDose.Data[[i]]$Time <- seq(from=10,by=10,to=400)


load("simdose.RData")

#############
# Simulation
#############

#                   sig11  S  OMEGA
SimDose.THETA <-  c(10 , 20 , .5)

SimDose.Data <- PSM.simulate(Model.SimDose, SimDose.Data, SimDose.THETA, dt=.1 , individuals=SimDose.Subj)

# View the Simulated datastructure
names(SimDose.Data[[1]])

# Plot of the simulations
par(mfcol=c(3,SimDose.Subj),mar = c(2, 4, 2, 2)+.1)
for(id in 1:SimDose.Subj) {
  plot(SimDose.Data[[id]]$Time , SimDose.Data[[id]]$Y,
         ylab="Observations", main=paste('individual',id,', eta:',round(SimDose.Data[[id]]$eta,3)))
  for(i in 1:2) {
    plot(SimDose.Data[[id]]$Time , SimDose.Data[[id]]$X[i,],type="l",
         ylab=paste('state',i))
    rug(SimDose.Data[[id]]$Time)
  }
  print(SimDose.Data[[id]]$eta)
}



###########
# Estimate
###########

parA <- list(LB=SimDose.THETA/50, Init=SimDose.THETA , UB=SimDose.THETA*50 )
fitA <- PSM.estimate(Model=Model.SimDose,Data=SimDose.Data,Par=parA,CI=T,trace=2)
fitA[1:3]

Model.SimDoseB <- Model.SimDose
Model.SimDoseB$ModelPar = function(THETA){
                        list(theta=list(sig11=0, S=THETA[1]), OMEGA=matrix(THETA[2]) ) }
parB <- list(LB=SimDose.THETA[2:3]/50, Init=SimDose.THETA[2:3] , UB=SimDose.THETA[2:3]*50 )

fitB <- PSM.estimate(Model=Model.SimDoseB,Data=SimDose.Data,Par=parB,CI=T,trace=2)
fitB[1:3]
#TEST for SDE
pchisq(2*(fitB$NegLogL-fitA$NegLogL), 3-2, lower.tail = FALSE) #p=0.005 - significant


#Plot individual LL with init param.
par(mfcol=c(3,2))
mp <- Model.SimDose$ModelPar(SimDose.THETA)
etavec <- seq(from=-2,to =2,length=50)
outvec <- vector(length=length(etavec))
gradvec <- vector(length=length(etavec))
for (id in 1:SimDose.Subj) {
  for (i in 1:length(etavec)) {
    outvec[i] <- IndividualLL.KF(etavec[i],mp$theta,mp$OMEGA,Model.SimDose,SimDose.Data[[id]])
    a <- IndividualLL.KF(etavec[i]+.001,mp$theta,mp$OMEGA,Model.SimDose,SimDose.Data[[id]])
    gradvec[i] <- (a-outvec[i])/.001
  }
  plot(etavec,outvec)
  abline(h=0)
}


###########
# Smoothing
###########

out <- PSM.smooth(Model = Model.SimDose, Data = SimDose.Data, THETA = fitA$THETA, subsample = 20)

# View the data structure
names(out[[1]])

#Plot of the smoothed estimates
par(mfcol=c(3,SimDose.Subj))
for(id in 1:SimDose.Subj) {
  plot(SimDose.Data[[id]]$Time , SimDose.Data[[id]]$Y,
         ylab="Observations", xlab=paste('individual',id))
  lines(out[[id]]$Time , out[[id]]$Xs[1,]/(V1*exp(out[[id]]$eta)))
  for(i in 1:2) {
    plot(out[[id]]$Time , out[[id]]$Xs[i,],type="l",
         ylab=paste('smooth state',i), xlab=paste('individual',id))
    lines(SimDose.Data[[id]]$Time , SimDose.Data[[id]]$X[i,],lty=2,col=4)
    rug(SimDose.Data[[id]]$Time)
    legend(400, y = max(out[[id]]$Xs[i,]), legend=c('simulation','smooth est'),lty=c(2,1),col=c(4,1),
           xjust=1,yjust=1,bty="n",cex=1,y.intersp=.7)
  }
}
