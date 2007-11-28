\name{sqrtm}
\alias{sqrtm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ <INTERNAL> Matrix square root }
\description{
  This is an internal function. Calculates the matrix square root using
  an eigenvalue decomposition. 
}
\usage{
sqrtm(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ A matrix }
}
}
\value{
  Returns the square root of the input matrix.
}
\author{ Stig B. Mortensen and Søren Klim }
\seealso{ \code{\link{PSM}}}
\examples{
(A <- matrix(c(5,1,1,8),nrow=2))
(sqrtA <- sqrtm(A))
sqrtA\%*\%sqrtA
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
