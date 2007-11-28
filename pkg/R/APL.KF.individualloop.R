`APL.KF.individualloop` <-
function (theta,OMEGA,Model,Data,GUIFlag=0) {
### NOTES -  requires: o$R, o$Yp, h(eta,theta)
  dimEta <- dim(OMEGA)[1]
  tmp <- dim(Data$Y)
  dimY <- tmp[1]
  dimN <- tmp[2]
  eGrad <- array(0,dim = c(dimY,dimEta,dimN))
  
  controllist <- list(trace=0,maxit=500,reltol=1e-7,ndeps=rep(1e-4,dimEta))
  
  if(GUIFlag>2) {
    controllist$trace <- 1  #higher number -> more detail
    controllist$REPORT <- 1 #report for every X iteration
  }
  
  out <- optim(par = rep(0,dimEta), fn = IndividualLL.KF , gr = IndividualLL.KF.gr,
               method = "BFGS", control = controllist, hessian = FALSE, 
               theta=theta, OMEGA=OMEGA, Model=Model, Data=Data)
  
  
  # Print optimization stats
  if(GUIFlag>2) {
    print(out$counts)
  }
  optimStat_i <- c(out$convergence, out$value, out$counts[1])
  
  phi <- Model$h(eta=out$par,theta=theta)
  o <- LinKalmanFilter(phi=phi, Model=Model, Data=Data, output=T)
  
  # Calculate stepsize in central diffenrence gradient
  stepSize <- 1E-5;
  
  # Create the e-grad
  for (p in 1:dimEta) {
    d <- rep(0,dimEta)
    d[p] <- stepSize
    
    # Forward difference
    phi.f <- Model$h(out$par+d,theta)
    eF <- Data$Y - LinKalmanFilter(phi=phi.f, Model=Model, Data=Data, output=T)$Yp

    # Backward difference
    phi.b <- Model$h(out$par-d,theta)
    eB <- Data$Y - LinKalmanFilter(phi=phi.b, Model=Model, Data=Data, output=T)$Yp;

    #Insert the calculated gradient.
    eGrad[,p,] <- (eF-eB)/(2*stepSize)
  }

  # Hessian approximation (19)
  h_Li <- matrix(0,dimEta,dimEta)
  for (q in 1:dimN) 
    h_Li = h_Li + t(C3(eGrad[,,q,drop=F]))%*%solve(o$R[,,q])%*%C3(eGrad[,,q,drop=F])
  h_Li <- - h_Li - solve(OMEGA);

  # RETURN neg. log. likelihood contribution
  list(   LiPart_i = .5*log(abs(det(h_Li))) + out$value,
             eta_i = out$par,
       optimStat_i = optimStat_i
       )
}

