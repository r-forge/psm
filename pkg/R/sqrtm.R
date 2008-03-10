`sqrtm` <-
function(A) {
  e <- eigen(A)
  return(e$vectors %*% diag(sqrt(e$values),nrow = length(e$values)) %*% t(e$vectors) )
}
