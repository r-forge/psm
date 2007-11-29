ModelCheck <- function(Model , Data , Par) {
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
  if( any(!( c("Matrices","X0","SIG","S","h","ModelPar","Dose")   %in% names(Model))) ) {
    print("Model doesn't contain the correct objects") ; return(FALSE) }
  # Model objects are functions

  for (name in c("Matrices","X0","SIG","S","h","ModelPar")) 
    if(!is.function(Model[[name]]))
      print(past("Model object",name,"is not a function"))

       
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
  if( any(!( c("Y","Time")  %in% names(Data))) ) {
    print("Data does not contain either Y og Time.")
    print("NOTE input to ModelCheck should be data for a single subject (eg. Data[[1]]).")
    return(FALSE) }


  
  # Calculate parameter phi
  Par <- Model$ModelPar(THETA=Par$Init)
  dimEta <- nrow(Par$OMEGA)
  phi <- Model$h(eta=rep(0,dimEta) , theta=Par$theta ,covar=Data$covar)

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
  if( any( eigen(SIG)$values <0 ) ) {
    print("Model$SIG is not positiv semidefinit") ; return(FALSE) }
  if( any( eigen(S)$values <0 ) ) {
    print("Model$S is not positiv semidefinit") ; return(FALSE) }

  
  
  # dimT
  dimT <- length(Data$Time) 
  if( dimT!=dim(Data$Y)[2]) {
    print("Data$Time and Data$Y doesn't match") ; return(FALSE) }
  if(ModelHasInput) {
    if(dimT!=dim(Data$U)[2]) {
      print("Data$Time and Data$U doesn't match") ; return(FALSE) }
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
  dimY <- nrow(matC)
  if (dimY != dim(Data$Y)[1] ){
    print("C and Data$Y doesn't match") ; return(FALSE) }
  if( dimY != nrow(S) | dimY != ncol(S) ) {
    print("S has incorrect size") ; return(FALSE) }

  if(ModelHasInput) {
    if (dimY != nrow(matD) ){
      print("C and D doesn't match") ; return(FALSE) }
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
    if( any(!(Model$Dose$Time %in% Data$Time))) {
      print("Dose times doesn't coincide with Data$Time") ; return(FALSE) }
      
    if( any( Model$Dose$State > dimX) ) {
      print("Dose states are larger than number of states") ; return(FALSE) }
  }

  return(TRUE)
}
