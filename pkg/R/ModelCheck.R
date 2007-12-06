ModelCheck <- function(Model , Data , Par, DataHasY=TRUE) {
  # Check the dimensions of the Model and Data
  #
  # Data is data for a single subject
  # Data$Time ; Data$Y ; Data$U


  if ( is.null(Data$U) ) {
    ModelHasInput <- FALSE
  } else if(any(is.na(Data$U))) {
    ModelHasInput <- FALSE
  } else {
    ModelHasInput <- TRUE
  }

  if (ModelHasInput) {
    Uk <- Data$U[1,]
  } else {
    Uk <- NULL
  }

  # Model contains the correct object
  if( any(!( c("Matrices","X0","SIG","S","h","ModelPar")   %in% names(Model))) ) {
    print("Model doesn't contain the correct objects") ; return(FALSE) }

  # Model objects are functions

  for (name in c("Matrices","X0","SIG","S","h","ModelPar")) 
    if(!is.function(Model[[name]])) {
      print(paste("Model object",name,"is not a function")); return(FALSE) }

       
#  if ( any(!(unlist(lapply(Model , function(x) {is.function(x)})))) ) {
#     print("Model objects are not all functions") ; return(FALSE) }

  # Check  Par
  if( !("Init" %in% names(Par)) ) {
    print("Init missing in Par") ; return(FALSE)
  }
  if( "LB" %in% names(Par) ) {
    #The parameter has bounds
    if( any( Par$LB > Par$Init) ) {
      print("Parameter LB is greater than Init") ; return(FALSE) }
    if( any( Par$UB < Par$Init) ) {
      print("Parameter UB is lower than Init") ; return(FALSE) }
  }

  # Check data
  if( any(!( "Time"  %in% names(Data))) ) {
    print("Individual data does not contain Time.") ; return(FALSE) }
  if(DataHasY) {
    if( any(!( "Y"  %in% names(Data))) ) {
      print("Individual data does not contain Y.") ; return(FALSE) }
    if( any(is.nan(Data$Y)) ) {
      print("Data$Y contains NaN. Use NA instead.");return(FALSE) }
  }
  if( any(c(is.nan(Data$Time),is.na(Data$Time))) ){
    print("Data$Time contains NA or NaN.");return(FALSE) }
  if (ModelHasInput && any(c(is.nan(Data$U),is.na(Data$U))) ){
    print("Data$U contains NA or NaN.");return(FALSE) }

  
  # Calculate parameter phi
  Parlist <- Model$ModelPar(THETA=Par$Init)
  if(!is.null(Parlist$OMEGA)) {
    dimEta <- nrow(Parlist$OMEGA)
    phi <- Model$h(eta=rep(0,dimEta) , theta=Parlist$theta ,covar=Data$covar)
  } else {
    phi <- Parlist$theta
  }

  # Matrices
  tmp   <- Model$Matrices(phi=phi)
  matA  <- tmp$matA
  matB  <- tmp$matB
  matC  <- tmp$matC
  matD  <- tmp$matD

  
  X0 <- Model$X0( Time=Data$Time[1], phi=phi, U=Uk)
  SIG <-  Model$SIG(  phi=phi)
  S <- Model$S(  phi=phi )

  # Test for positive semidefinit
  if(dim(SIG)[1]!=dim(SIG)[2] || dim(SIG)[1]==0) {
    print("Model$SIG is not a square matrix") ; return(FALSE) }
  if(dim(S)[1]!=dim(S)[2] || dim(S)[1]==0) {
    print("Model$S is not a square matrix") ; return(FALSE) }
      
  if( any( eigen(SIG)$values <0 ) ) {
    print("Model$SIG is not positiv semidefinit") ; return(FALSE) }
  if( any( eigen(S)$values <0 ) ) {
    print("Model$S is not positiv semidefinit") ; return(FALSE) }

  
  # dimT
  dimT <- length(Data$Time)
  if(DataHasY) {
    if( dimT!=dim(Data$Y)[2]) {
      print("Data$Time and Data$Y does not have same length") ; return(FALSE) }
  }
  if(ModelHasInput) {
    if(dimT!=dim(Data$U)[2]) {
      print("Data$Time and Data$U does not have same length") ; return(FALSE) }
    } # ModelHasInput

  # dimX
  dimX <- nrow(matA)
  if( dimX != ncol(matA) ) {
    print("A is not square [dimX dimX]") ; return(FALSE) }
  if( dimX!=ncol(matC) ) {
    print("A and C dimensions doesn't match") ; return(FALSE) }
  if( dimX!=length(X0) ) {
    print("X0 has incorrect size") ; return(FALSE) }
  if( dimX != nrow(SIG) | dimX != ncol(SIG) ) {
    print("SIG has incorrect size") ; return(FALSE) }
  
  if(ModelHasInput) {
    if( dimX!=nrow(matB) ) {
      print("A and B dimensions doesn't match") ; return(FALSE) }
    } # ModelHasInput

  # dimY
  if(DataHasY) {
    dimY <- nrow(matC)
    if (dimY != dim(Data$Y)[1] ){
      print("C and Data$Y doesn't match") ; return(FALSE) }
    if( dimY != nrow(S) | dimY != ncol(S) ) {
      print("S has incorrect size") ; return(FALSE) }
  }

  if(ModelHasInput && DataHasY) {
    if (dimY != nrow(matD) ){
      print("Y and D doesn't match") ; return(FALSE) }
  }
  
  # dimU
  if( ModelHasInput) {
    dimU <- nrow(Data$U)
    if( dimU != ncol(matB)) {
      print("Data$U and B doesn't match") ; return(FALSE) }
    if( dimU != ncol(matD)) {
      print("Data$U and D doesn't match") ; return(FALSE) }
  }
  
  # Dose
  if( "Dose" %in% names(Model) ) {
    MD <- Model$Dose
    if( any(!(MD$Time %in% Data$Time))) {
      print("Dose times doesn't coincide with Data$Time") ; return(FALSE) }
      
    if( any( MD$State > dimX) ) {
      print("Dose states are larger than number of states") ; return(FALSE) }

    if( !all( (c(length(MD$State),length(MD$Amount))- length(MD$Time))==0  ) ) {
      print("Dose: Elements Time, State and Amount not of same length") ; return(FALSE)}
  }

  return(TRUE)
}
