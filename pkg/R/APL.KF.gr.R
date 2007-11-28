`APL.KF.gr` <-
function (THETA,Model,Pop.Data,LB=NULL,UB=NULL,GradSTEP=1e-4,GUIFlag=0) {
  # Forward gradient function for APL.KF

  
  L <- length(THETA)
  APL <- rep(0,length=(L+1))
  GRAD <- rep(0,length=L)
  
  TP <- matrix(0,nrow=L+1,ncol=L)
  for (i in 1:L) {
    TP[i, ] <- THETA
    TP[i,i] <- THETA[i] + GradSTEP * ( abs(THETA[i])  + GradSTEP );
  }
  TP[L+1,] <- THETA
  
  for ( i in 1:(L+1) )
    APL[i] <- APL.KF(TP[i,],Model=Model,Pop.Data=Pop.Data,LB=LB,UB=UB,GUIFlag=0)
        
  # Calculate Gradient and insert
  for (i in 1:L)
    GRAD[i] = (APL[i]-APL[L+1])/( TP[i,i]-THETA[i] ) 

  if(GUIFlag>1) {
    if(!is.null(LB))
      print(c(" THETA  =", signif(invlogit(THETA,LB,UB),5)),q=F)
    print(c("<THETA> =", signif(THETA,5)),q=F)
    print(c("<GR>    =", signif(GRAD,5)),q=F)
  }
   
  GRAD
}

