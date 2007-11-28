`PSM.simulate` <- 
function(Model,THETA,dt,Ulist=NULL,Tlist,individuals=1) {
  # LinSimulation
  # FUNCTION FOR SIMULATION OF LINEAR MIXED EFFECTS MODEL USING SDE'S.
  #--------
  #INPUT
  # Model - Model list
  # THETA - Model parameters
  # dt    - time step in discretization
  # Ulist - list of input for each individual, each element may be
  #           NA: no input assumed (for any individual)
  #           dim(U[[i]])=(?,1): constant input assumed
  #           dim(U[[i]])=(?,n): where n=[t(end)-t(start)]/dt+1
  # Tlist - List of sample times for each individual
  # individuals - Number of individuals to simulate
  #--------
  #OUTPUT
  # List containing Xlist, Ylist, Tlist, Ulist & eta


  Result <-  vector(mode="list",length=individuals)
  # Xlist <- Ylist <- vector(mode="list",length=individuals)
  OMEGA <- Model$ModelPar(THETA)$OMEGA
  theta <- Model$ModelPar(THETA)$theta

  if(!is.null(OMEGA)) {
    dimEta <- dim(OMEGA)[1]
    eta <- sqrtm(OMEGA) %*% matrix(rnorm(individuals*dimEta),nrow=dimEta,ncol=individuals)
  } else {
    eta <- NULL
  }


  for (i in 1:individuals) {
    print(paste("Simulating individual: ",i))
    if(!is.null(OMEGA)) {
      phi <- Model$h(eta[,i],theta)
    }else {
      phi <- theta
    }

    SampleTime <- Tlist[[i]]
    len <- length(SampleTime)
    t <- seq(from = SampleTime[1],to=SampleTime[len],by=dt)
    TimeSpan <- SampleTime[len]-SampleTime[1]
    n <- TimeSpan/dt + 1
    
    # find where SampleTime are equal to t
    where <- abs(t(matrix(rep(SampleTime,n),nrow=len,ncol=n))-t) < 1000*.Machine$double.eps
    
    
    if(sum(where)!=len)
      stop(simpleError("All sample times must belong to t(0)+n*dt, where n is an integer."))
    


    # Check for INPUT and subsample U
    if( is.null(Ulist) ) { #check if U exists.
      ModelHasInput <- FALSE
    }
    else if ( any(is.na(Ulist)) ) { #check if it contains any NA
      ModelHasInput <- FALSE
    }
    else {
      ModelHasInput <- TRUE
    }

    if( ModelHasInput) {
      RepIdx <- rep( 1:length(SampleTime) , times=c(diff(SampleTime)/dt,1))
      U <- Ulist[[i]][ , RepIdx ]
      Ustart <- U[,1,drop=F]
    } else {
      Ustart <- NA
    }
    
    ModelHasDose <- "Dose" %in% names(Model)
    
    # Create Matrices
    tmpM  <- Model$Matrices(phi=phi)
    matA  <- tmpM$matA
    matB  <- tmpM$matB
    matC  <- tmpM$matC
    matD  <- tmpM$matD 


    # Initial States
    InitialState <- Model$X0(phi=phi)
    SIG <- Model$SIG(phi=phi)  

    # Dimensions
    dimX <- nrow(InitialState)
    dimY <- nrow(matC)
    
    # Is A singular ???
    rankA <- qr(matA)$rank
    singA <- (rankA<dimX)
  
    X <- array(0,dim=c(dimX,n))
    Y <- array(0,dim=c(dimY,n))
  
  # Variables in the for loop
    DIAGDIMX  <- diag(1,dimX)
    DIAGRANKA <- diag(1,rankA)

  #Gaussian noise
    eW <-   matrix(rnorm(n*dimX),nrow=dimX,ncol=n)
    eObs <- matrix(rnorm(n*dimX),nrow=dimY,ncol=n)
    
    X[,1] <- InitialState
    
    
    for (k in 1:n) {
 
      if(ModelHasInput) 
        Uk <- U[,k,drop=F]

      # observation
      if(ModelHasInput) {
        Y[,k] <- matC %*% X[,k,drop=F] + matD%*%Uk + sqrtm(Model$S(phi=phi)) %*% eObs[,k]
      } else 
      Y[,k] <- matC %*% X[,k,drop=F] + sqrtm(Model$S(phi=phi)) %*% eObs[,k]

      # Add dose after measurement is taken at Time[k]
      if(ModelHasDose) {
        idxD = which(t[k]==Model$Dose$Time)
        if(length(idxD)==1) {
          X[Model$Dose$State[idxD],k] <- X[Model$Dose$State[idxD],k] + Model$Dose$Amount[idxD]
        }
      }
      
 
    # Abort state pred if finished
      if (k == n) break

    # State prediction
      tmp <- dt * rbind(cbind(-matA , SIG%*%t.default(SIG)) ,
                        cbind( matrix(0,nrow=dimX,ncol=dimX) , t.default(matA) ))
      tmp <- mexp(tmp)
      PHI <- t.default( tmp[(dimX+1):(2*dimX),(dimX+1):(2*dimX),drop=F])
      IntExpASIG <- PHI %*% tmp[1:dimX,(dimX+1):(dimX*2),drop=F]

    # Different formulaes depending on A and INPUTS
      if( !singA ) {
      # Special case #3: Non-singular A, zero order hold on inputs.
        X[,k+1] <- { if( ModelHasInput) {
          PHI%*%X[,k,drop=F] + solve.default(matA)%*%(PHI-DIAGDIMX)%*%matB%*%Uk
        } else {
          PHI%*%X[,k,drop=F] } }
      } else {
      # Special case #1: Singular A, zero order hold on inputs.
        if(!ModelHasInput) {
          X[,k+1] <- PHI %*% X[,k,drop=F]
        }
        else  {
          if(rankA==0) { #only the zero matrix has rank 0, thus A=0: dx=b*u*dt
            X[,k+1] <- X[,k] +  matB*Uk*dt
          } else {
          
                  # matA is singular and has INPUT and dimX>1
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
            IntExpAtildeS[(rankA+1):dimX,(rankA+1):dimX] <- diag(1,dimX-rankA)*dt

                  # Insert State prediction CTSM  (1.60)
            X[,k+1] <- PHI %*% X[,k]+ Ua %*% IntExpAtildeS %*% t.default(Ua) %*% matB %*% Uk
          } # end if (A is zero matrix) else ...
        } # end else
      } # end else
      
    # Wiener noise
      ew <- SIG %*% (sqrt(dt) * eW[,k+1])
      X[,k+1] <- X[,k+1] + ew
    
    } #end for

    # SubSampling
    idx <- which(where %*% rep(T,len) == 1)

    Result[[i]] <- list(X=matrix(X[,idx],nrow=dimX),Y=matrix(Y[,idx],nrow=dimY),
                        Time=SampleTime,U=Ulist[[i]],eta=eta)

  } #end individual loop
  
  # list(Xlist=Xlist,Ylist=Ylist,Tlist=SampleTlist,Ulist=Ulist,eta=eta)
  return(Result)
}



#matrix square root
sqrtm <- function(A) {
  e <- eigen(A)
  e$vectors %*% diag(sqrt(e$values),nrow = length(e$values)) %*% t(e$vectors)
}



