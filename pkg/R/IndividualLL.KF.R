`IndividualLL.KF` <-
function (eta,theta,OMEGA,Model,Data){
### NOTES -  requires: o$negLogLike, o$Yp, h(eta,theta)

  phi <- Model$h(eta,theta)

  # run the KF one time, to evaluate negative log-likelihood
  negLogLike <- LinKalmanFilter( phi=phi, Model=Model , Data=Data )
  eta <- matrix(eta,ncol=1)

  #Return a posteriori negative log likelihood
  negLogLike + 0.5*t(eta)%*%solve(OMEGA)%*%eta + 0.5*log( det(2*pi*OMEGA))

}

