PSM.template <- function(Linear=FALSE,dimX=2,dimY=3,dimU=4,dimEta=5,Dose=TRUE,file="") {
  
  str <- "MyModel <- vector(mode=\"list\")"
  
  if(Linear) {
    str <- paste(str,"MyModel$Matrices=function(phi) {",sep="\n")
    str <- paste(str,"  list(",sep="\n")
    str <- paste(str,"\n       matA=matrix(c(  ), nrow=",dimX,", ncol=",dimX,"),",sep="")
    if(dimU)
    str <- paste(str,"\n       matB=matrix(c(  ), nrow=",dimX,", ncol=",dimU,"),",sep="")
    str <- paste(str,"\n       matC=matrix(c(  ), nrow=",dimY,", ncol=",dimX,"),",sep="")
    if(dimU)
    str <- paste(str,"\n       matD=matrix(c(  ), nrow=",dimX,", ncol=",dimU,")",sep="")
    str <- paste(str,"       )",sep="\n")
    str <- paste(str,"}",sep="\n")
  } else {
    str <- paste(str,"MyModel$Functions <- ",sep="\n")
    str <- paste(str,"  list(",sep="\n")
    str <- paste(str,"       f = function(x,u,time,phi) {",sep="\n")
    str <- paste(str,"\n         matrix(c(  ), nrow=",dimX,", ncol=1)",sep="")  
    str <- paste(str,"       },",sep="\n")
    str <- paste(str,"       df = function(x,u,time,phi) {",sep="\n")
    str <- paste(str,"\n         matrix(c(  ), nrow=",dimX,", ncol=",dimX,")",sep="")  
    str <- paste(str,"       },",sep="\n")
    str <- paste(str,"       g = function(x,u,time,phi) {",sep="\n")
    str <- paste(str,"\n         matrix(c(  ), nrow=",dimY,", ncol=1)",sep="")  
    str <- paste(str,"       },",sep="\n")
    str <- paste(str,"       dg = function(x,u,time,phi) {",sep="\n")
    str <- paste(str,"\n         matrix(c(  ), nrow=",dimY,", ncol=",dimX,")",sep="")  
    str <- paste(str,"       }",sep="\n")
    str <- paste(str,"       )",sep="\n")
  }

  str <- paste(str,"MyModel$h = function(eta,theta,covar) {",sep="\n")
  str <- paste(str,"  phi <- theta",sep="\n")
  str <- paste(str,"  phi",sep="\n")
  str <- paste(str,"}",sep="\n")

  if(Linear) {
    str <- paste(str,"MyModel$S = function(phi) {",sep="\n")
  } else {
    str <- paste(str,"MyModel$S = function(u,time,phi) {",sep="\n")
  }
  str <- paste(str,"\n  matrix(c(  ), nrow=",dimY,", ncol=",dimY,")",sep="")
  str <- paste(str,"}",sep="\n")
   
  if(Linear) {
    str <- paste(str,"MyModel$SIG = function(phi) {",sep="\n")
  } else {
    str <- paste(str,"MyModel$SIG = function(u,time,phi) {",sep="\n")
  }
  str <- paste(str,"\n  matrix(c(  ), nrow=",dimX,", ncol=",dimX,")",sep="")
  str <- paste(str,"}",sep="\n")
  
  str <- paste(str,"MyModel$X0 = function(Time,phi,U) {",sep="\n")
  str <- paste(str,"\n  matrix(c(  ), nrow=",dimX,", ncol=1)",sep="")
  str <- paste(str,"}",sep="\n")
  
  str <- paste(str,"MyModel$ModelPar = function(THETA) {",sep="\n")
  str <- paste(str,"  list(theta=list(  )",sep="\n")
  if(dimEta)
    str <- paste(str,",\n       OMEGA=matrix(c(  ), nrow=",dimEta,", ncol=",dimEta,")\n       ",sep="")
  str <- paste(str,")\n}",sep="")

  if(Dose)
    str <- paste(str,"MyModel$Dose = list(Time=c(  ), State=c(  ), Amount=c(  ))",sep="\n")

  str <- paste("\n",str,"\n\n",sep="\n")

  cat(str,file=file)
  
}
