`invlogit` <-
function(y,xmin,xmax)
  exp(y)/(1+exp(y))*(xmax-xmin)+xmin

