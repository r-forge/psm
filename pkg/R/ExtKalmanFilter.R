`ExtKalmanFilter` <-
function( phi, Model, Data, outputInternals=FALSE) {

  # Extract Data components
  Time  <- Data[["Time"]] 
  Y     <- Data[["Y"]]

  # Check for INPUT and set it
  if(is.null(Data[["U"]])) { #check if U exists.
    ModelHasInput <- FALSE
  }
  else if ( any(is.na(Data[["U"]])) ) { #check if it contains any NA
    ModelHasInput <- FALSE
  }
  else {
    ModelHasInput <- TRUE
  }
  U     <- if( !ModelHasInput) { NA } else { Data[["U"]] }

  # Set Uk
  if(ModelHasInput) {
    Uk <- U[,1,drop=FALSE]
  } else {Uk <- NA}

  ModelHasDose <- "Dose" %in% names(Model)

  f  <- Model$Functions$f
  df <- Model$Functions$df
  g  <- Model$Functions$g
  dg <- Model$Functions$dg

  # Calculate Init state and initial SIG
  InitialState  <- Model$X0(Time=Time[1], phi=phi, U = Uk)
  SIG           <- Model$SIG(u=Uk,time=Time[1],phi=phi)
  Ax0 <- df(x=InitialState ,u=Uk,time=Time[1],phi=phi)  ######### Er DETTE RIGTIGT??


  dimT  <- length(Time)     # Time is vector -> use length
  dimY  <- nrow(Y)          # Dimensionality of observations
  dimU  <- ifelse(ModelHasInput,nrow(U),0)
  dimX  <- nrow(InitialState)

  # P0 CTSM MathGuide page 19 (1.118) and page 8 (1.49)
  PS  <- 1.0;
  tau <- Time[2] - Time[1]

  # Dimensions Pintegral:[2*dimX 2*dimX]
  tmp <-  rbind( cbind(-Ax0 , SIG%*%t.default(SIG)) ,
                cbind( matrix(rep(0,2*dimX),nrow=dimX,ncol=dimX) , t.default(Ax0) ))*tau
  
  # Use Matrix package to compute Matrix exponential
  Pint  <- matexp(tmp)
  
  PHI0   <- t.default(Pint[(dimX+1):(2*dimX),(dimX+1):(2*dimX),drop=FALSE])
  P0    <- PS*PHI0 %*% Pint[1:dimX,(dimX+1):(2*dimX),drop=FALSE]

  #----------------------------------------------------------
  # Init matrices used in the Kalman filtering
  #----------------------------------------------------------
  Xp      <- array(NA,c(dimX,dimT))
  Xf      <- array(NA,c(dimX,dimT))
  Yp      <- array(NA,c(dimY,dimT))
  KfGain  <- array(NA,c(dimX,dimY,dimT))
  Pf      <- array(NA,c(dimX,dimX,dimT))
  Pp      <- array(NA,c(dimX,dimX,dimT))
  R       <- array(NA,c(dimY,dimY,dimT))

  # Insert initial estimates into matrices
  Pp[,,1] <- P0
  Xp[,1]  <- InitialState


  ######################
  # Loop over timepoints
  ######################

  negLogLike <- 0
 
  for(k in 1:dimT) {

    # Set Uk
    if(ModelHasInput) {
      Uk <- U[,k,drop=FALSE]
    } else {Uk <- NA}                    
    
    
    ######################
    # Update
    ######################
    
    # Y_Hat, Prediction
    Yp[,k] <- g(x=Xp[,k,drop=FALSE],u=Uk,time=Time[k],phi=phi)
    
    # Find the h-derivative in this point
    # dg(X,U,t,phi)
    C <- dg(x=Xp[,k,drop=FALSE],u=Uk,time=Time[k],phi=phi)
    S <- Model$S(u=Uk,time=Time[k],phi=phi)

    # Uncertainty on Measurement.
    R[,,k] <- C%*%Pp[,,k,drop=FALSE]%*%t(C) + S
    
    # Kalman gain
    KfGain[,,k] <- Pp[,,k,drop=FALSE]%*%t(C)%*%solve(R[,,k,drop=FALSE])
    
    # Innovation
    e <- Y[,k,drop=FALSE] - Yp[,k,drop=FALSE]
    
    # Updating Equations
    Xf[,k] <- Xp[,k,drop=FALSE] + KfGain[,,k,drop=FALSE]%*%e
    
    # Pf(:,:,k) = Pp(:,:,k) - KfGain(:,:,k)*R(:,:,k)*KfGain(:,:,k)';
    Pf[,,k] <- CutThirdDim(Pp[,,k,drop=FALSE]) -
      CutThirdDim(KfGain[,,k,drop=FALSE])%*%C%*%CutThirdDim(Pp[,,k,drop=FALSE])

    # Add contribution to negLogLike
    negLogLike = negLogLike + 0.5 * ( log(det(2*pi*CutThirdDim(R[,,k,drop=FALSE])))
      + t(e)%*%solve(R[,,k])%*%e )

    # Abort state pred if finished
    if(k==dimT) break

    
    ######################
    # Prediction
    ######################
    
    # Create Z combined variable
    # Upper triangle of P
    Index <- NULL
    for (p in 1:dimX) {
      Index <- c(Index , (1:p) + dimX*(p-1) )
    }
    tmpP <- Pf[,,k,drop=FALSE]
    tmpP <- tmpP[Index]
    Z <- c( Xf[,k] , tmpP)

    dSystemPred <- function(t,y,parms) {
      # Evaluate dX
      X  <- matrix(y[1:dimX],ncol=1)
      dX <- f(x=X,u=Uk,time=t,phi=phi)
      # Evaluate dP
      tmpP <- matrix(0,ncol=dimX,nrow=dimX)
      tmpP[Index] <- y[-(1:dimX)]
      if(dimX>1) {
        tmpP <- tmpP + t(tmpP) - diag(diag(tmpP))
      }
      Ax <- df(x=X,u=Uk,time=t,phi=phi)
      SIGx <- Model$SIG(u=Uk,time=t,phi=phi)
      ### MathGuide (1.79)
      dP <- Ax%*%tmpP + tmpP*t(Ax) + SIGx%*%t(SIGx) 
      dP <- dP[Index]
      # Return dX and dP
      list(c( dX, dP))
    }

    # Prediction of Z
    timevec <- c(Time[k], Time[k+1]) #!should be subsampled for smoothing!!  
    ZOUT <- lsoda(y=Z, times=timevec, func=dSystemPred, parms=NULL, rtol=1e-6, atol=1e-6)
    
#    %save foreward state estimates for smoothing.
#    if(nargout==2)
#        sub_foreward(i).time = dummyT;
#        sub_foreward(i).state = ZOUT(:,1:dimX);
#    end
    
    # convert back to X,Pk
    Xp[,k+1] <- ZOUT[length(timevec),1+(1:dimX)] #first col is time
    
    tmpP <- matrix(0,ncol=dimX,nrow=dimX)
    tmpP[Index] <- ZOUT[length(timevec),-(1:(dimX+1))]
    if(dimX>1) {
      tmpP <- tmpP + t(tmpP) - diag(diag(tmpP))
    }
    Pp[,,k+1] <- tmpP
    
  } #end loop over observations

  if(outputInternals) {
    return( list( negLogLike=negLogLike,Time=Time,Xp=Xp, Xf=Xf, Yp=Yp, KfGain=KfGain,Pf=Pf, Pp=Pp, R=R))
  } else {
    return(as.vector(negLogLike))
  }

}
