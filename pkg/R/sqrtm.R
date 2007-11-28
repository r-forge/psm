`sqrtm` <-
function(A) {
  e <- eigen(A)
  e$vectors %*% diag(sqrt(e$values),nrow = length(e$values)) %*% t(e$vectors)
}
