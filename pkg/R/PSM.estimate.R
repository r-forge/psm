`PSM.estimate` <-
function(Model,Data,Par,CI=FALSE,trace=0,optimizer="optim", controllist=NULL,fast=TRUE,...) {

  dimS <- length(Data)
  for(i in 1:dimS) {
    check <- ModelCheck(Model,Data[[i]],Par)
    if(!check$ok) {
      print(paste("Error occured using data for individual",i))
      break
    }
  }
  if(!check$ok) stop(paste("Input did not pass model check."))
  Linear <- check$Linear

  if(trace>0)
    cat(ifelse( Linear, "* Linear model *\n", "* Non-linear model *\n"))
  
  # Check for fast option and Singular A
  if(Linear && fast) {
    #Get the ModelParameters
    tmp       <- Model$ModelPar(Par$Init)
    if( is.null(tmp$OMEGA) ) {
       # OMEGA IS NULL, but covariates can be present      
        tmpPhi <- Model$h( eta=NULL , theta=tmp$theta , covar=Data[[1]]$covar)                
    } else {
      # OMEGA is not null
      tmpDimEta <- dim(tmp$OMEGA)[1] 
      tmpPhi    <- Model$h( eta=rep(0, tmpDimEta) , theta=tmp$theta , covar=Data[[1]]$covar)
    }
    tmpMat    <- Model$Matrices(tmpPhi)
    matA  <- tmpMat$matA
    
    # Check for Model INPUT (U)
    if(is.null(Data[[1]][["U"]])) { #check if U exists.
      ModelHasInput <- FALSE
    } else if ( any(is.na(Data[[1]][["U"]])) ) { #check if it contains any NA
      ModelHasInput <- FALSE
    } else {
      ModelHasInput <- TRUE
    }
    tmpU     <- if( !ModelHasInput) { NA } else { Data[[1]][["U"]] }

    # Set Uk
    if(ModelHasInput) {
      Uk <- tmpU[,1,drop=FALSE]
    } else {Uk <- NA}
    
    tmpdimX <- nrow( Model$X0(Time=Data[[1]]$Time[1], phi=tmpPhi, U = Uk) )
    
    rankA <- qr(matA)$rank
    singA <- (rankA<tmpdimX)
    
    if( singA) {
      cat("Unable to use option \"fast\" on singular A matrix - Switching to fast=FALSE \n") 
      fast=FALSE 
    }
  }

  T0 <- proc.time()[3]

  if(!is.null(Par$LB)) {
    if(trace>1) cat("Using logit transformation of parameters \n")
    Par$Init <- logit(Par$Init,Par$LB,Par$UB) }

  # controllist
  if( is.null(controllist)) {
    # The user did not supply a controllist for the optimizer
    controllist <- list(maxit=100, trace=trace, REPORT=1 )
    # put parameters on similar scale if bounds are missing
    if(is.null(Par$LB))
      controllist$parscale <- abs(Par$Init)+1e-3
  }
  
  switch( EXPR = optimizer,
    optim={
    if(trace>1)  cat( "Using Optimizer: \t optim\n")
    out <- optim(par = Par$Init, fn = APL.KF ,
                 gr = APL.KF.gr, method = "BFGS",
                 control = controllist, hessian = CI,
                 Model=Model, Pop.Data=Data, LB=Par$LB, UB=Par$UB,
                 GUIFlag=trace,fast=fast,Linear=Linear,...)
    NegLogL     <- out$value
    ParEstimate <- out$par
    if(CI) Hess <- out$hessian
    
    },
    nlm={ 
      if(trace>1)  cat( "Using Optimizer: \t nlm\n")

      # typsize=Par$Init,stepmax=(.1*abs(Par$Init))+1e-3,          
      out <- nlm(f=APL.KF, p=Par$Init, hessian=CI, print.level=trace,
                 Model=Model, Pop.Data=Data, LB=Par$LB,
                 UB=Par$UB, GUIFlag=trace,fast=fast,Linear=Linear,...)
                 
        NegLogL     <- out$minimum 
        ParEstimate <- out$estimate
        if(CI) Hess <- out$hessian
    },
    "Otherwise" = {       stop( "Optimizer not recognized")     } 
    )
    
    
  
  ci <- CI
  if(CI) { 
    ci <- matrix(c(
                   out$par-1.96*sqrt(diag(solve(Hess))),
                   out$par,
                   out$par+1.96*sqrt(diag(solve(Hess)))),nrow=3,byrow=TRUE)
    rownames(ci) <- c("Lower CI95","MLE","Upper CI95")
    colnames(ci) <- names(Par$Init)
  } 
  
  if(!is.null(Par$LB)) {
    ParEstimate = invlogit(ParEstimate,Par$LB,Par$UB)
    if(CI)
      for (i in 1:3) {
        ci[i,] <- invlogit(ci[i,],Par$LB,Par$UB)
      }
  } 
  
  totaltime = proc.time()[3]-T0
  
  if(trace>0) {
    minutes = floor(totaltime/60)
    tid <- paste("Runtime:  ",minutes,":",round(totaltime-60*minutes,2),"",sep="")
    print(tid,quote=FALSE)
  }

  list(NegLogL = NegLogL, THETA = ParEstimate, CI = ci, opt = out, sec = totaltime)

}

