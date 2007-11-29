`LinKalmanFilter` <-
function( phi , Model , Data , echo=F, outputInternals=FALSE) {
# Linear Continous Kalman Filter for a State Space formulation
#     Input: 
#     Model           type: list
#       $Matrices     function  input:  phi,u
#       $X0           function  input:  phi,u,t
#       $SIG          function  input:  phi,u,t
#
#     Data            type: list
#       $Time         matrix [1 dimT]
#       $Y            matrix [dimY dimT]
#       $U            matrix [dimU dimT]
#
#     phi             Parameter vector
#     echo            TRUE or FALSE   Display informaiton during execution
#
#     Problem definition
#     dx = (Ax + Bu)dt + SIG(phi,u,t) db
#     y = Cx + Du + e             e ~ N(0,S(phi,u,t))
#
#     Developer Notes:
#
#     Important dimensions
#       dimT    number of timepoints (observations)
#       dimX    number of states
#       dimU    dimension of input
#       dimY    output dimension
#
#     A:[dimX dimX] ; B:[dimX dimU] ; SIG:[dimX dimX]
#     C:[dimY dimX] ; D:[dimY dimU] ; S:[dimY,dimY]
#

    # Print current value of phi .. Handy in minimizations.
  if(echo) cat("phi:" , paste(round(as.double(phi),2) , "\t"))

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
      Uk <- U[,1,drop=F]
    } else {Uk <- NA}

    ModelHasDose <- "Dose" %in% names(Model)
     
    tmp   <- Model$Matrices(phi=phi)
    matA  <- tmp$matA
    matB  <- tmp$matB
    matC  <- tmp$matC
    matD  <- tmp$matD 
  
    # Calculate Init state and initial SIG
    InitialState  <- Model$X0(Time=Time[1], phi=phi, U = Uk)
    SIG           <- Model$SIG(phi=phi)


    dimT  <- length(Time)     # Time is vector -> use length
    dimY  <- nrow(Y)          # Dimensionality of observations
    dimU  <- ifelse(ModelHasInput,nrow(U),0)
    dimX  <- nrow(InitialState)

    # Is A singular ???
    rankA <- qr(matA)$rank
    singA <- (rankA<dimX)
    
    # P0 CTSM MathGuide page 19 (1.118) and page 8 (1.49)
    PS  <- 1.0;
    tau <- Time[2] - Time[1]

    # Dimensions Pintegral:[2*dimX 2*dimX]
    tmp <-  rbind( cbind(-matA , SIG%*%t.default(SIG)) ,
                  cbind( matrix(rep(0,2*dimX),nrow=dimX,ncol=dimX) , t.default(matA) ))*tau
#print(tmp)
    # Use Matrix package to compute Matrix exponential
    Pint  <- mexp(tmp)

    PHI0   <- t.default(Pint[(dimX+1):(2*dimX),(dimX+1):(2*dimX),drop=F])
    P0    <- PS*PHI0 %*% Pint[1:dimX,(dimX+1):(2*dimX),drop=F]

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


    #----------------------------------------------------------
    # Kalman Filtering for Linear Time-Invariant Model
    #----------------------------------------------------------

    # Init the Negative LogLikelihood
    negLogLike <- 0

    # Variables in the for loop
    DIAGDIMY  <- diag(1,dimY)
    DIAGDIMX  <- diag(1,dimX)
    DIAGRANKA <- diag(1,rankA)

    # Loop over timepoints
    for(k in 1:dimT) {
        # Does the observation has missing values?
        ObsIndex  <- which(!is.na(Y[,k]))
        E         <- DIAGDIMY[ObsIndex,,drop=F]

        # Set Uk
        if(ModelHasInput) {
          Uk <- U[,k,drop=F]
        } else {Uk <- NA}                

        # Create output prediction
        Yp[ObsIndex,k] <- {
                  if(ModelHasInput) {
                  E %*% (matC%*%Xp[,k,drop=F] + matD%*%Uk)
                  } else { 
                  E%*%matC%*%Xp[,k,drop=F]} }

        # Output prediction covariance    
        S     <- Model$S(phi=phi)
        R[ObsIndex,ObsIndex,k]  <- E%*%matC%*%Pp[,,k]%*%t.default(matC)%*%t.default(E) + E%*%S%*%t.default(E)

  
        if(length(ObsIndex)>0) { #there must be at least one obs for updating.
          # Kalman gain
          KfGain[,ObsIndex,k] <- Pp[,,k]%*%t.default(matC) %*% t.default(E) %*% solve.default(R[ObsIndex,ObsIndex,k])

          # Updating
          e       <- Y[ObsIndex,k,drop=F]-Yp[ObsIndex,k,drop=F]
          KFg     <- CutThirdDim(KfGain[,ObsIndex,k,drop=F])
          Xf[,k]  <- Xp[,k,drop=F] + KFg%*%e
          Pf[,,k] <- Pp[,,k] - KFg %*% R[ObsIndex,ObsIndex,k] %*% t.default(KFg)
  
          # Add contribution to negLogLike: CTSM page 3. (1.14)
          tmpR        <-  CutThirdDim(R[ObsIndex,ObsIndex,k,drop=F])
          tmp        <-  determinant.matrix(tmpR)
     
          negLogLike  <- negLogLike + .5*( log(tmp$sign*exp(tmp$modulus)) + t.default(e)%*%solve.default(tmpR)%*%e)
        } else { #no observations, update not available.
          Xf[,k] <- Xp[,k]
          Pf[,,k] <- Pp[,,k]
        }

        if(ModelHasDose) {
          idxD = which(Time[k]==Model$Dose$Time)
          if(length(idxD)==1) {
            Xf[Model$Dose$State[idxD],k] <- Xf[Model$Dose$State[idxD],k] + Model$Dose$Amount[idxD]
          }
        }
        
        # Abort if negLogLike is Inf or -Inf
        if(is.infinite(negLogLike)) {
          if(echo) cat("\t -LL: " ,  round(negLogLike,3), "\n")
          return(negLogLike)
          }
        
        # Abort state pred if finished
        if(k==dimT) break
        

        # State prediction
        tau   <- Time[k+1]-Time[k]

        # Try to optimize
        tmp   <- tau* rbind( cbind(-matA , SIG%*%t.default(SIG)) ,
                  cbind( matrix(0,nrow=dimX,ncol=dimX) , t.default(matA) ))
        
        tmp   <- mexp(tmp)

        # CTSM (1.48)
        PHI   <- t.default(tmp[(dimX+1):(2*dimX),(dimX+1):(2*dimX),drop=F])
        # CTSM (1.49)
        IntExpASIG <- PHI %*% tmp[1:dimX,(dimX+1):(dimX*2),drop=F]
        # CTSM (1.45)
        Pp[,,k+1] <- PHI %*% Pf[,,k] %*% t.default(PHI) + IntExpASIG

        # Different formulaes depending on A and INPUTS
        if( !singA ) {
            # Special case #3: Non-singular A, zero order hold on inputs.
            Xp[,k+1] <- { if( ModelHasInput) {
                              PHI%*%Xf[,k,drop=F]+ solve.default(matA)%*%(PHI-DIAGDIMX)%*%matB%*%Uk
                              } else {
                              PHI%*%Xf[,k,drop=F] } }
        } else {

              # Special case #1: Singular A, zero order hold on inputs.
              if(!ModelHasInput) {
                  Xp[,k+1] <- PHI %*% Xf[,k]
                }
              else  {
                  # matA is singular and has INPUT
                  # CTSM Mathguide Special Case no.1 page 10
                  Ua <- svd(matA)$u

                  PHITilde      <- t.default(Ua) %*% PHI %*% Ua
                  PHITilde1     <- PHITilde[1:rankA,1:rankA,drop=F]
                  # PHITilde2     <- PHITilde[1:rankA,(rankA+1):dimX,drop=F]
                  # PHITilde1Inv  <- solve(PHITilde1)

                  ATilde    <- t.default(Ua) %*% matA %*% Ua
                  ATilde1   <- ATilde[1:rankA,1:rankA,drop=F]
                  ATilde2   <- ATilde[1:rankA,(rankA+1):dimX,drop=F]
                  ATilde1Inv <- solve.default(ATilde1)

                  IntExpAtildeS <- matrix(NA,dimX,dimX)
                  # Insert upper left part of matrix [1:rankA 1:rankA]
                  IntExpAtildeS[1:rankA , 1:rankA] <- ATilde1Inv %*% (PHITilde1-DIAGRANKA)
                  # Lower left part
                  IntExpAtildeS[(rankA+1):dimX , 1:rankA]   <- 0
                  # Upper Right
                  IntExpAtildeS[1:rankA,(rankA+1):dimX] <- ATilde1Inv %*%
                      (IntExpAtildeS[1:rankA,1:rankA]-DIAGRANKA*tau)%*%ATilde2
                  # Lower right
                  IntExpAtildeS[(rankA+1):dimX,(rankA+1):dimX] <- diag(1,dimX-rankA)*tau

                  # Insert State prediction CTSM  (1.60)
                  Xp[,k+1] <- PHI %*% Xf[,k]+ Ua %*% IntExpAtildeS %*% t.default(Ua) %*% matB %*% Uk

                  } # end else
            }    # end else
        } #end for

        # Complete negLogLike  with second part CTSM page 3 (1.14)
        #negLogLike  <- negLogLike + .5*dimT*dimY*log(2*pi)
        negLogLike  <- negLogLike + .5*( dimT*dimY - sum(is.na(Y)) )*log(2*pi)

        if(echo) cat("\t -LL: " ,  round(negLogLike,0) , "\n")
        if(outputInternals) {
          return( list( negLogLike=negLogLike,Xp=Xp, Xf=Xf, Yp=Yp, KfGain=KfGain,Pf=Pf, Pp=Pp, R=R))
        } else {
          return(as.vector(negLogLike)) }
}





