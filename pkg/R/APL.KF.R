`APL.KF` <-
function (THETA,Model,Pop.Data,LB=NULL,UB=NULL,GUIFlag=0,longOutput=F) {
### NOTES -  requires: Model$ModelPar(), 
  
  if (!is.null(LB)) {
    THETA <- invlogit(THETA,LB,UB)
  }

  OMEGA <- Model$ModelPar(THETA)$OMEGA
  theta <- Model$ModelPar(THETA)$theta
  

  # Dimensions
  dimS <- length(Pop.Data)
  ifelse(is.null(OMEGA),dimEta <- 1, dimEta <- dim(OMEGA)[1])

  # Init
  etaList <- matrix(0,nrow=dimEta,ncol=dimS)
  optimStat <- matrix(0,nrow=3,ncol=dimS)
  LiPart <- matrix(0,nrow=1,ncol=dimS)
  
  if(GUIFlag>1) 
    starttime <- proc.time()[3]
  

  for (i in 1:dimS) {
    if(GUIFlag>2)
      print(paste('Individual', i))
    if(!is.null(OMEGA)) {
      result <- APL.KF.individualloop(theta=theta,OMEGA=OMEGA,Model=Model,Data=Pop.Data[[i]],GUIFlag=GUIFlag) 
      LiPart[i] <- result$LiPart_i
      etaList[,i] <- result$eta_i
      optimStat[,i] <- result$optimStat_i
    }
    else {
      LiPart[i] <- LinKalmanFilter( phi=theta, Model=Model , Data=Pop.Data[[i]] )
      etaList[,i] <- NaN
      optimStat[,i] <- NaN
    }
      
  }
  
  if(GUIFlag>1) {
    totaltime = proc.time()[3]-starttime
    minutes = floor(totaltime/60)
    tid <- paste("  (",minutes,":",round(totaltime-60*minutes,2),")",sep="")
    print(c(" -logL  =", signif(sum(LiPart),10),tid),q=F)

  }
  
  if(longOutput) {
    list(negLogLike=sum(LiPart),etaList=etaList,optimStat=optimStat)
  } else {
    #The return variable - neg. Log. Likelihood
    sum(LiPart)
  }
}

