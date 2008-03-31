`APL.KF.individualloop` <-
function (theta,OMEGA,Model,Data,GUIFlag=0,fast=TRUE) {
### NOTES -  requires: o$R, o$Yp, h(eta,theta)
  dimEta <- dim(OMEGA)[1]
  tmp <- dim(Data$Y)
  dimY <- tmp[1]
  dimN <- tmp[2]
  eGrad <- array(0,dim = c(dimY,dimEta,dimN))
  
  controllist <- list(trace=0,maxit=500)
  
  if(GUIFlag>2) {
    controllist$trace <- GUIFlag-2  #higher number -> more detail
    controllist$REPORT <- 1 #report for every X iteration
  }
  if(0) { #unconstrained
    controllist$reltol <- 1e-7
    out <- optim(par = rep(0,dimEta), fn = IndividualLL.KF , gr = IndividualLL.KF.gr,
                 method = "BFGS", control = controllist, hessian = FALSE, 
                 theta=theta, OMEGA=OMEGA, Model=Model, Data=Data, fast=fast)
  } else { #constrained
    controllist$factr <- 1e8
    out <- optim(par = rep(0,dimEta), fn = IndividualLL.KF, gr = IndividualLL.KF.gr,
                 method = "L-BFGS-B", lower = -4*sqrt(diag(OMEGA)),upper = 4*sqrt(diag(OMEGA)),
                 control = controllist, hessian = FALSE, 
                 theta=theta, OMEGA=OMEGA, Model=Model, Data=Data, fast=fast)
  }
  
  
  # Print optimization stats
  if(GUIFlag>2) {
    print(out$counts)
  }
  optimStat_i <- c(out$convergence, out$value, out$counts[1])
  
  phi <- Model$h(eta=out$par,theta=theta,covar=Data$covar)
  o <- LinKalmanFilter(phi=phi, Model=Model, Data=Data, output=TRUE, fast=fast)
  
  # Calculate stepsize in central diffenrence gradient
  stepSize <- 1E-5;
  
  # Create the e-grad
  for (p in 1:dimEta) {
    d <- rep(0,dimEta)
    d[p] <- stepSize
    
    # Forward difference
    phi.f <- Model$h(out$par+d,theta,covar=Data$covar)
    eF <- Data$Y - LinKalmanFilter(phi=phi.f, Model=Model, Data=Data, outputInternals=TRUE, fast=fast)$Yp

    # Backward difference
    phi.b <- Model$h(out$par-d,theta,covar=Data$covar)
    eB <- Data$Y - LinKalmanFilter(phi=phi.b, Model=Model, Data=Data, outputInternals=TRUE, fast=fast)$Yp;

    #Insert the calculated gradient.
    eGrad[,p,] <- (eF-eB)/(2*stepSize)
  }

  # Hessian approximation (19)
  h_Li <- matrix(0,dimEta,dimEta)
  for (q in 1:dimN) {
    ObsIndex  <- which(!is.na(Data$Y[,q]))
    if(length(ObsIndex)>0) {
      eG <- CutThirdDim(eGrad[ObsIndex,,q,drop=FALSE])
      h_Li <- h_Li + t.default(eG) %*%  solve(o$R[ObsIndex,ObsIndex,q]) %*% eG
    }
  }
  h_Li <- - h_Li - solve(OMEGA);

  # RETURN neg. log. likelihood contribution
  list(   LiPart_i = .5*log(abs(det(h_Li))) + out$value,
             eta_i = out$par,
       optimStat_i = optimStat_i
       )
}

