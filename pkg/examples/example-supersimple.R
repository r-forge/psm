library(PSM)


mod1 <- vector(mode="list")
mod1$Matrices=function(phi) {
  list(
       matA=matrix(c(-1*phi$k),ncol=1),
       matC=matrix(c(1/phi$V),nrow=1)
       )
}
mod1$X0 = function(Time,phi,U) {
  A0 <- phi$A0
  matrix(c(A0),ncol=1)
}
mod1$S = function(phi) {
  matrix(phi[["S"]])
}
mod1$ModelPar = function(THETA) {
  S <- 20
  A0 <- 1500
  list(theta=list(k=THETA['k'],V=THETA['V'], A0 = A0, S = S),
       OMEGA=matrix(.2))
}
mod1$SIG = function(phi) {
  matrix(10)
}
mod1$h = function(eta,theta,covar) {
  phi = theta
  phi$V = theta$V*exp(eta[1])
  phi
}

names(mod1)
TimeVec <- c(0,1,2,3,5,7,10,13,16,20,25)
PrepData = list( list(Time = TimeVec),list(Time = TimeVec) )

MyTHETA = c(k = 0.08, V = 15)

SimData <- PSM.simulate(mod1,PrepData,MyTHETA,deltaTime = .1)


names(SimData[[1]])
plot(SimData[[1]]$Time,SimData[[1]]$Y)
lines(SimData[[1]]$longTime,SimData[[1]]$longX/MyTHETA['V']/exp(SimData[[1]]$eta))
points(SimData[[2]]$Time,SimData[[2]]$Y,col=2)
lines(SimData[[2]]$longTime,SimData[[2]]$longX/MyTHETA['V']/exp(SimData[[2]]$eta),col=2)


par <- list(LB = MyTHETA/20,
            Init = c(k=.005,V=100),
            UB = MyTHETA*20)

fit <- PSM.estimate(mod1,SimData,par,CI=TRUE)
fit[1:3]

sm <- PSM.smooth(mod1,SimData,fit$THETA,subs=10)
names(sm[[1]])
lines(sm[[1]]$Time,sm[[1]]$Ys,col=3)

##### Penalty function.


mini = .1
maxi = 10
xvec = seq(mini,maxi,length=1000)
lambda = 1e-4
dx = 1e-30
dx=0
p = lambda*(mini/(xvec-mini+dx) + maxi/(maxi-xvec+dx))
range(p)
plot(xvec,p,type="l",ylim=c(0,.001))



##########################################################
### Non-linear model


mod2 <- mod1[c("X0","h","ModelPar")]
mod2$Functions <- 
  list(
       f = function(x,u,time,phi) {
         -1*phi$k*x
       },
       df = function(x,u,time,phi) {
         -1*phi$k
       },
       g = function(x,u,time,phi) {
         x/phi$V
       },
       dg = function(x,u,time,phi) {
         1/phi$V
       }
       )
mod2$S = function(u,time,phi) {
  matrix(phi[["S"]])
}
mod2$SIG = function(u,time,phi) {
  matrix(10)
}
MyPar <- mod2$ModelPar(MyTHETA)
myphi <- mod2$h(0,MyPar$theta,NA)
unlist(myphi)

source(file="~/PSM/PSM/R/ExtKalmanFilter.R")
source(file="~/PSM/PSM/R/PSM.estimate.R")
source(file="~/PSM/PSM/R/APL.KF.R")
source(file="~/PSM/PSM/R/APL.KF.gr.R")
source(file="~/PSM/PSM/R/IndividualLL.KF.R")
source(file="~/PSM/PSM/R/APL.KF.individualloop.R")
source(file="~/PSM/PSM/R/IndividualLL.KF.gr.R")
source(file="~/PSM/PSM/R/PSM.smooth.R")


ExtKalmanFilter( myphi, mod2, SimData[[1]] )
LinKalmanFilter( myphi, mod1, SimData[[1]] )


mod1.no.omega <- mod1
mod2.no.omega <- mod2
mod1.no.omega$ModelPar = function(THETA) {
  list(theta=list(k=THETA['k'],V=THETA['V'], A0 = 1500, S = 20), OMEGA=NULL)
}
mod2.no.omega$ModelPar <- mod1.no.omega$ModelPar

# Linear, no-OMEGA
fit1.no.omega <- PSM.estimate(mod1.no.omega,SimData,par)
fit1.no.omega[1:2]

# Non-Linear, no-OMEGA
fit2.no.omega <- PSM.estimate(mod2.no.omega,SimData,par)
fit2.no.omega[1:2]

# Linear, with-OMEGA
fit1 <- PSM.estimate(mod1,SimData,par)
fit1[1:2]

# Non-Linear, with-OMEGA
fit2 <- PSM.estimate(mod2,SimData,par,trace=1)
fit2[1:2]

