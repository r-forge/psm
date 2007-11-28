
k1 = 0.05; k2 = 0.05; ke = 0.07;
Model.SimDose <- list(
                      Matrices = function(phi) {
                        kei <- phi[["kei"]]
                        matA <- matrix( c(-(k1+kei) ,  k2 ,   
                                                 k1 , -k2 ) ,nrow=2,byrow=T)
                        matC <- matrix(c(1,0),nrow=1)
                        list(matA=matA,matC=matC)
                      },
                      X0 = function(phi) {
                        matrix(0,nrow=2)
                      },
                      SIG = function(phi) { 
                        diag( c(0,0) )
                      },
                      S = function(phi) {
                        matrix(phi[["S"]])
                      },
                      h = function(eta,theta) {
                        phi <- theta
                        phi[["kei"]] <- theta[["ke"]]*exp(eta[1])
                        phi
                      },
                      ModelPar = function(THETA){
                        list(theta=list(ke=THETA[1], S=THETA[2]),
                             OMEGA=diag(c(1)) )
                      },
                      Dose = list(
                        Time = c(30,150),
                        State = c(1, 1),
                        Amount = c(1000,1500)
                        )
                      )

# Create Simulation Timeline 
SimDose.Subj <- 1
SimDose.Time <- vector(mode="list",length=NoOfSubjects)
for (i in 1:NoOfSubjects) 
  SimDose.Time[[i]] <- seq(from=15,by=15,length=30)

SimDose.Time


detach(package:PSM)
library(PSM,lib.loc="~/PSM/Rpackages/gridterm")

#                   ke   S  OMEGA
SimDose.THETA <-  c(0.03 , 0 , 0)
Model.SimDose$Dose$Time = c(30,180)
Model.SimDose$Dose$Amount = c(1500,1500)
SimDose.Data <- PSM.simulate(Model.SimDose, SimDose.THETA, dt=.1 , Tlist=SimDose.Time ,individuals=SimDose.Subj)


par(mfcol=c(2,2))
for(id in 1:SimDose.Subj) {
  for(i in 1:2) {
    plot(SimDose.Data[[id]]$Time , SimDose.Data[[id]]$X[i,],type="l",
         ylab=paste('state',i), xlab=paste('individual',id))
    rug(SimDose.Data[[id]]$Time)
  }
}

