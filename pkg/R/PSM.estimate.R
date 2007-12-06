`PSM.estimate` <-
function(Model,Data,Par,CI=F,trace=0,optimizer="optim", controllist=NULL) {

  ok <- TRUE
  dimS <- length(Data)
  for(i in 1:dimS) {
    ok <- ModelCheck(Model,Data[[i]],Par)
    if(!ok) {
      print(paste("Error occured using data for individual",i))
      break
    }
  }
  if(!ok) stop(paste("Input did not pass model check."))
    
  T0 <- proc.time()[3]

  if(!is.null(Par$LB)) {
    Par$Init <- logit(Par$Init,Par$LB,Par$UB) }

  # controllist
  if( is.null(controllist)) {
    # The user did not supply a controllist for the optimizer
    controllist <- list(maxit=100, abstol=1e-5, trace=trace, REPORT=1 )
    # put parameters on similar scale if bounds are missing
    if(is.null(Par$LB))
      controllist$parscale <- abs(Par$Init)+1e-3
  }
  
  if(trace>1)  cat( "Using Optimizer:", optimizer , "\n")

  if(optimizer=="optim") {
    out <- optim(par = Par$Init, fn = APL.KF ,
                 gr = APL.KF.gr, method = "BFGS",
                 control = controllist, hessian = CI,
                 Model=Model, Pop.Data=Data, LB=Par$LB, UB=Par$UB,
                 GUIFlag=trace)
    NegLogL <- out$value
    ParEstimate <- out$par
    if(CI) Hess <- out$hessian
    
  } else if(optimizer=="nlm") {
    StepLength <- 0.5
    if(is.null(Par$LB)) {
      out <- nlm(f=APL.KF, p=Par$Init, hessian=CI, print.level=trace,
                 typsize=Par$Init,
                 stepmax=(StepLength*abs(Par$Init)) +1e-3,
                 Model=Model, Pop.Data=Data, LB=Par$LB,
                 UB=Par$UB, GUIFlag=trace)
    } else {
      cat( "Optimizer not recognized -> Using Optimizer: nlm \n")
      out <- nlm(f=APL.KF, p=Par$Init, hessian=CI, print.level=trace,
                 Model=Model, Pop.Data=Data, LB=Par$LB,
                 UB=Par$UB, GUIFlag=trace)
    }
    
    NegLogL <-out$minimum 
    ParEstimate <- out$estimate
    if(CI) Hess <- out$hessian
  }
  
  ci <- CI
  if(CI) { 
    ci <- matrix(c(
                   out$par-1.96*sqrt(diag(solve(Hess))),
                   out$par,
                   out$par+1.96*sqrt(diag(solve(Hess)))),nrow=3,byrow=T)
    rownames(ci) <- c("Lower CI95","MLE","Upper CI95")
  } 
  
  if(!is.null(Par$LB)) {
    THETA = invlogit(out$par,Par$LB,Par$UB)
    if(CI)
      for (i in 1:3) {
        ci[i,] <- invlogit(ci[i,],Par$LB,Par$UB)
      }
  } else {
    THETA = out$par
  }
  totaltime = proc.time()[3]-T0
  
  if(trace>0) {
    minutes = floor(totaltime/60)
    tid <- paste("Runtime:  ",minutes,":",round(totaltime-60*minutes,2),"",sep="")
    print(tid,q=F)
  }

  list(NegLogL = out$value, THETA = THETA, CI = ci, opt = out, sec = totaltime)

}

