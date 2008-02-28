`invlogit` <-
function(y,xmin=0,xmax=1)
  exp(y)/(1+exp(y))*(xmax-xmin)+xmin

