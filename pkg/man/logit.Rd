\name{(inv)logit}
\alias{logit}
\alias{invlogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ <INTERNAL> (inv)logit transformation.} 
\description{
  These are internal functions. Gives logit and invlogit
  transformation of a vector. 
}
\usage{
logit(x, xmin, xmax)
invlogit(y, xmin, xmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{ Vectors on normal and logit scale }
  \item{xmin, xmax}{ Vectors of upper and lower bounds on normal scale}
}

\value{
  Returns a vector of transformed values.
}

\author{ Stig B. Mortensen and Søren Klim }
\seealso{ \code{\link{PSM}}}
%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
