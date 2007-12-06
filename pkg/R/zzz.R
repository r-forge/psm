.First.lib <- function(lib, pkg) {

  library.dynam("PSM", pkg, lib)
  
  loaded <- suppressWarnings(require(mexp))
    if(!loaded) {
     	print("Searching for package in lib.loc")
    	require(mexp,lib.loc=lib)
    	}
    
}
