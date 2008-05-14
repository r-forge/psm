`IndividualLL.KF` <-
function (eta,theta,OMEGA,Model,Data,fast=TRUE,Linear){
### NOTES -  requires: o$negLogLike, o$Yp, h(eta,theta)
  
  phi <- Model$h(eta,theta,covar=Data$covar)

  # run the KF one time, to evaluate negative log-likelihood
  if(Linear) {
    negLogLike <- LinKalmanFilter( phi=phi, Model=Model , Data=Data , fast=fast)
  } else {
    negLogLike <- ExtKalmanFilter( phi=phi, Model=Model , Data=Data )
  } 
  eta <- matrix(eta,ncol=1)

  #Return a posteriori negative log likelihood
  negLogLike + 0.5*t(eta)%*%solve(OMEGA)%*%eta + 0.5*log( det(2*pi*OMEGA))

}

