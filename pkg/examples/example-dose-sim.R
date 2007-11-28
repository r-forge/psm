
detach(package:PSM)
library(PSM,lib.loc="~/PSM/Rpackages/gridterm")


k1 = 0.05; k2 = 0.05; ke = 0.07;
Model.SimDose <- list(
                      Matrices = function(phi, U) {
                        kei <- phi[["kei"]]
                        matA <- matrix( c(-(k1+kei) ,  k2 ,   
                                                 k1 , -k2 ) ,nrow=2,byrow=T)
                        matC <- matrix(c(1,0),nrow=1)
                        list(matA=matA,matC=matC)
                      },
                      X0 = function(Time, phi, U) {   ############## remove Time, U
                        matrix(0,nrow=2)
                      },
                      SIG = function(Time,phi, U) {   ############## remove Time, U
                        diag( c(0,0) )
                      },
                      S = function(Time , phi = phi, U ) {   ############## remove Time, U
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
                      }
                      )

# Create Simulation Timeline 
SimDose.Subj <- 2
SimDose.Time <- vector(mode="list",length=NoOfSubjects)
for (i in 1:NoOfSubjects) 
  SimDose.Time[[i]] <- seq(from=15,by=15,length=30)

SimDose.Time

SimDose.THETA <-  c(0.07 , 100 , 1)

Sim.Data <- PSM.simulate(Model.SimDose, SimDose.THETA, dt=.1 , Tlist=SimDose.Time ,individuals=SimDose.Subj)
