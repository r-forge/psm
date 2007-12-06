\name{CutThirdDim}
\alias{CutThirdDim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ <INTERNAL> Remove third dimension of dim-attribute. }
\description{
  This is an internal function. Removes third and higher dimensions of
  dim-attribute for an array and thus creating a matrix.
}
\usage{
CutThirdDim(a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{ An array with length(dim(a))>=2}
}

\value{
  Returns a matrix.
}
\author{ Stig B. Mortensen and Søren Klim }
\seealso{ \code{\link{PSM}}}
\examples{
tmp  <- array(1:4,dim=c(2,2,1))
tmp2 <- matrix(2:5,nrow=2,ncol=2)
err  <-try(tmp\%*\%tmp2,silent=TRUE)
err[1]  
CutThirdDim(tmp)\%*\%tmp2
print("No error!")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
