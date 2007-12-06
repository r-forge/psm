# Test Case
# Originating from CTSM
# Linear Time Invarient - Heat Case

rm(list=ls())

detach(package:PSM)

library(PSM,lib.loc="~/PSM/Rpackages/gridterm")

# Load the Data and Variables
tmpData <- read.table("Heat_Data.csv",sep=";", col.names=c("TIME","Te","Ti","Q"))

Time=tmpData$TIME
Y=t(matrix(tmpData[,c("Q")]))
U=t(as.matrix(tmpData[,c("Te","Ti")]))
Pop.Data <- list( list(Time=tmpData$TIME,
                Y=t(matrix(tmpData[,c("Q")])),
                U=t(as.matrix(tmpData[,c("Te","Ti")]))) )

HeatModel <- list(
                  Matrices = function(phi=NA) {
                    G1  <- phi[["G1"]] ; G2  <- phi[["G2"]]
                    H1  <- phi[["H1"]] ; H2  <- phi[["H2"]] ; H3  <- phi[["H3"]]
                    tmp <- list(
                                matA = matrix( c(-1*(1/H1+1/H2)/G1,1/(G1*H2),1/(G2*H2) , -1*(1/H2+1/H3)/G2 ) , ncol=2, byrow=T),
                                matB = diag( c(1/(G1*H1) , 1/(G2*H3) ) ),
                                matC = matrix( c(0,-1/H3) ,nrow=1),
                                matD = matrix( c(0,1/H3)  ,nrow=1))
                    return(tmp)
                  },
                  X0 = function(Time=NA,phi=NA,U=NA) {
                    tmp    <- phi[["X01"]]
                    tmp[2] <- phi[["X02"]]
                    return(matrix(tmp,ncol=1) )} ,
                  SIG = function(phi=NA) {
                    return( diag( c(phi[["SIG11"]],phi[["SIG22"]])))} ,
                  S = function(phi=NA) {
                    return( matrix(phi[["S"]])) } ,
                  h = function(eta,theta,covar=NULL) {
                    phi <- theta
                    return(phi) } ,
                  ModelPar = function(THETA){
                    return(list(theta=list( G1=THETA[1],G2=THETA[2],
                                  H1=THETA[3],H2=THETA[4],H3=THETA[5],
                                  SIG11=THETA[6], SIG22=THETA[7], S=THETA[8],
                                  X01=THETA[9], X02=THETA[10]),
                                OMEGA=NULL))},
                  Dose=NULL
                  )

names(HeatModel)

# Parameter estimation
# Initial guess from CTSM
# THETA OBJ             G1,   G2,  H1,  H2,  H3, SIG11,SIG22,    S, X01,  X02
# CTSM starting guess fails in this implementation
par1 <- list(LB   = c(   0,    0,   0,   0,   0,     0,    0,    0,   10,   20),
             Init = c( 100,   50,   1,   2,  .5,   .01,  .01,  .01,   15,   25),
             UB   = c( 200,  100,   2,   5,   1,     1,    1,    1,   20,   30)
             )

par1$Init <-        c( 100 ,  50,   1,   2,  .5,  .001, .001, .001,   13,   25)
par1$UB <- par1$LB <- NULL


# Check the Model
ModelCheck( Model=HeatModel , Data=Pop.Data[[1]], Par=par1)


# -------------------------------------------------------------
# Test Linear Kalman Filter with CTSM estimated parameters
# -------------------------------------------------------------

# CTSM returns -LL= -623 
CTSMphi <- c( 1.3134E+01,2.5330E+01,1.0394E+02,9.6509E-01,2.0215E+00,4.9320E+01,5.0929E-01,7.3779E-08,2.6951E-09,1.0330E-02)
names(CTSMphi) <- c("X01","X02","G1","H1","H2","G2","H3","SIG11","SIG22","S")
Ob1 <- LinKalmanFilter( phi=CTSMphi , Model=HeatModel , Data=Pop.Data[[1]] , echo=T, outputInternals=TRUE)

# Validation plot versus Data
D <- Pop.Data[[1]]
plot(D$Time , D$Y,col="red")
points( D$Time , Ob1$Yp, pch="+",col="blue")

# -------------------------------------------------------------
# Minimizers
# -------------------------------------------------------------
                                        # Test Run of initial parameters
phi <- par1
names(phi$Init) <- c("G1","G2","H1","H2","H3","SIG11","SIG22","S","X01","X02")
phi$Init

# Perform minimization with 2 different optimizers
Rprof()
Min1 <- PSM.estimate(Model=HeatModel,Data=Pop.Data,Par=par1,CI=TRUE,trace=2,optimizer="nlm")
Rprof(NULL)
summaryRprof()

Min2 <- PSM.estimate(Model=HeatModel,Data=Pop.Data,Par=par1,CI=TRUE,trace=2,optimizer="optim")

cat( "nlm: "   , Min1$sec, "\t ", Min1$opt$minimum, "\n")
cat( "optim: " , Min2$sec, "\t ", Min2$opt$value, "\n")


# -------------------------------------------------------------
# Smoother 
# -------------------------------------------------------------

SmoothObj <- PSM.smooth(Model=HeatModel, Data=Pop.Data, THETA=CTSMphi, subsample=0,trace=1)


D <- SmoothObj[[1]]
names(D)

Idx <- 200:500
plot( D$Time[Idx], D$Xs[1,Idx] , type="n" )
for(i in 1:2) polygon( c(D$Time[Idx],rev(D$Time[Idx])) , c(D$Xs[i,Idx],rev(D$Xs[i,Idx]))+sqrt( abs(c(D$Ps[i,i,Idx], - rev(D$Ps[i,i,Idx])))),col=4)
for(i in 1:2) lines( D$Time[Idx], D$Xs[i,Idx], type="l",lwd=2)



# Load the predicted Data from CTSM
CTSM.Val.Data <- read.table("CTSM_Heat_Xp.csv",sep=";", col.names=c("Time","Xp1","Xp2","SDX1","SDX2","Y1","SDY1"))

names(D)


# Measurements Validation
plot(Pop.Data[[1]]$Time, Pop.Data[[1]]$Y  )
lines( CTSM.Val.Data$Time, CTSM.Val.Data$Y1, col="red")
Idx <- !is.na(D$Yp)
lines( D$Time[Idx] , D$Yp[Idx] , col="blue")

# State validation
plot(D$Time[-(1:10)], D$Xp[1,-(1:10)],type="l",col="red")
lines(D$Time[-(1:10)], D$Xp[2,-(1:10)],type="l",col="red")

plot(CTSM.Val.Data$Time , CTSM.Val.Data$Xp2,type="l",col="blue",ylim=range(CTSM.Val.Data[,c("Xp1","Xp2")]))
lines(CTSM.Val.Data$Time , CTSM.Val.Data$Xp1,type="l",col="blue")


